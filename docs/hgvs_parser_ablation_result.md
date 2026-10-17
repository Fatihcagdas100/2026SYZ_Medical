# HGVS ayrıştırıcı denetimi ve özellik ablasyonu

Tarih: 18 Ekim 2026

## Amaç ve veri politikası

Bu adım iki soruyu yanıtlar:

1. Ham ClinVar `Name` alanından türetilen HGVS sınıfları sözdizimsel olarak
   beklenen örneklerle ve kaydedilmiş özellik tablosuyla tutarlı mı?
2. Gene-disjoint validation başarımındaki artış hangi özellik basamağında
   oluşuyor?

Bütün model karşılaştırmaları aynı 257.035 train ve 55.068 validation kaydı,
aynı lojistik regresyon ayarları, `seed=2027` ve sabit `0.5` karar eşiğiyle
yapıldı. 55.093 satırlık test bölümü geliştirme veri çerçevesine yüklenmedi ve
bu analizde test metriği üretilmedi.

## Ayrıştırıcı denetimi

Gene-disjoint train/validation içindeki 312.103 kayıt yeniden ayrıştırıldı.
Kaydedilmiş özellik tablosuyla karşılaştırılan altı alanın hiçbirinde uyuşmazlık
bulunmadı:

- `hgvs_sequence_type`
- `dna_change_type`
- `coding_region`
- `splice_offset_category`
- `protein_effect`
- `has_protein_annotation`

Missense, synonymous, frameshift, stop-gained/changed, protein anotasyonsuz ve
belirsiz/absent sınıflarının her birinden sabit SHA-256 sıralamasıyla 10'ar ham
`Name` örneği seçildi. Toplam 60 örneğin tamamında üretim ayrıştırıcısından ayrı
basit işaret kontrolü beklenen sözdizimini buldu.

Bu sonuç, kodun aynı girdiden aynı özellikleri ürettiğini ve seçilen örneklerde
temel HGVS işaretlerini doğru yakaladığını gösterir. Klinik anlam doğruluğunu,
referans transkript seçimini veya HGVS standardının bütün uç durumlarını
doğrulamaz; bunlar için uzman incelemesi gerekir.

Yerel inceleme dosyası:
`reports/generated/hgvs_parser_audit/sampled_raw_name_examples.tsv`.

## Üç basamaklı özellik ablasyonu

Ana karşılaştırma modeli açık VCF eksiklik bayrağını kullanmayan lojistik
regresyondur.

| Özellik profili | ROC-AUC | AP | Duyarlılık | Özgüllük | F1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Yalnız yapısal | 0.7630 | 0.7516 | 0.6147 | 0.7889 | 0.6678 |
| Yapısal + HGVS biyolojik sonuç | 0.9825 | 0.9810 | 0.8912 | 0.9561 | 0.9196 |
| Biyolojik sonuç + açık gösterim alanları | 0.9824 | 0.9810 | 0.8902 | 0.9568 | 0.9194 |

Yapısal profile HGVS biyolojik sonuçlarının eklenmesi:

- ROC-AUC: `+0.219460`
- AP: `+0.229418`
- dengeli doğruluk: `+0.221834`
- F1: `+0.251795`

Açık `hgvs_sequence_type` ve `has_protein_annotation` alanlarının son basamakta
eklenmesi:

- ROC-AUC: `-0.000123`
- AP: `-0.000063`
- dengeli doğruluk: `-0.000175`
- F1: `-0.000227`

Son basamağın pratik katkısı gözlenmemiştir. Ancak orta profilde
`protein_effect=no_protein_annotation` kategorisi bulunur; dolayısıyla anotasyon
yokluğu dolaylı olarak hâlâ temsil edilebilir. Son basamak bütün
anotasyon-varlığı sinyalini değil, açık gösterim alanlarının ek katkısını ölçer.

## Tam modelin validation alt grupları

| Alt grup | n | Benign | Pathogenic | Precision | Recall | ROC-AUC | AP | Özgüllük |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Tüm varyantlar | 55.068 | 28.495 | 26.573 | 0.9505 | 0.8902 | 0.9824 | 0.9810 | 0.9568 |
| Missense | 7.075 | 3.884 | 3.191 | 0.5262 | 0.3149 | 0.5696 | 0.5078 | 0.7670 |
| Synonymous | 5.338 | 5.292 | 46 | 0.0000 | 0.0000 | 0.5311 | 0.0114 | 1.0000 |
| Frameshift | 11.430 | 69 | 11.361 | 0.9940 | 1.0000 | 0.5582 | 0.9949 | 0.0000 |
| Stop-gained/changed | 7.773 | 24 | 7.749 | 0.9969 | 1.0000 | 0.4667 | 0.9967 | 0.0000 |
| Protein anotasyonu var | 32.203 | 9.444 | 22.759 | 0.9491 | 0.8994 | 0.9658 | 0.9843 | 0.8838 |
| Protein anotasyonu yok | 22.865 | 19.051 | 3.814 | 0.9596 | 0.8351 | 0.9746 | 0.9411 | 0.9930 |

## Yorum

Genel skor artışı bütün protein sonuç sınıflarında güçlü ayrım öğrenildiği
anlamına gelmez. Missense grubunda ROC-AUC yalnızca `0.5696` ve duyarlılık
`%31,49`dur; sıradaki modelleme adımının ana hedefi budur.

Synonymous grubunda yalnızca 46 pozitif, frameshift grubunda yalnızca 69 negatif,
stop grubunda yalnızca 24 negatif vardır. Bu aşırı sınıf dengesizliği nedeniyle
AP, duyarlılık ve özgüllük tek başına güçlü ayrım kanıtı değildir. Özellikle
frameshift ve stop gruplarında model bütün veya hemen hemen bütün kayıtları
patojenik tahmin etmektedir.

Sonuçlar HGVS biyolojik sonuç sınıflarının genel gene-disjoint validation
performansına büyük katkı verdiğini, fakat missense ayrımının çözülmediğini
gösterir. Bu bir validation bulgusudur; klinik geçerlilik iddiası değildir.

## Yeniden çalıştırma

```powershell
.\.venv\Scripts\python.exe scripts\audit_hgvs_parser.py
.\.venv\Scripts\python.exe scripts\run_hgvs_ablation.py
```

Otomatik testler:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```
