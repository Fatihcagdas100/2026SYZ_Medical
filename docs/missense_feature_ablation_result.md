# Missense biyolojik özellik ablasyonu

Tarih: 18 Ekim 2026

## Amaç

Önceki analizde genel HGVS modeli güçlü görünmesine rağmen 7.075 missense
validation kaydında ROC-AUC `0.5696`, recall `%31,49` bulunmuştu. Bu deney,
gen veya protein kimliği kullanmadan yalnızca aminoasit değişiminin biyolojik
niteliğinin bu gruptaki ayrımı geliştirip geliştirmediğini ölçer.

## Üretilen özellikler

Standart `p.Arg41Gly` benzeri missense gösteriminden şu alanlar çıkarıldı:

- referans ve alternatif aminoasit,
- yönlü aminoasit değişim çifti (`R>G` gibi),
- Grantham uzaklığı,
- BLOSUM62 ikame skoru,
- Kyte–Doolittle hidropati farkı ve mutlak farkı,
- temel kalıntı sınıfı ve yük sınıfı geçişleri.

Protein konumu, gen sembolü, transkript kimliği ve ham `Name` model tablosuna
alınmadı. Özellikler hedef etikete bakılmadan üretildi. Toplam 47.708 missense
kaydın tamamı standart aminoasit çifti olarak ayrıştırıldı.

Denetim sırasında `Xaa` ile bilinmeyen alternatif aminoasit belirten üç train
kaydının eski genel kalıp tarafından missense sınıfına alındığı bulundu. Genel
HGVS ayrıştırıcısı yalnızca 20 standart aminoasidi missense kabul edecek şekilde
düzeltildi; bu üç kayıt `uncertain_or_absent` sınıfına taşındı. İlgili özellik
tabloları ve geliştirme raporları yeni hash'lerle yeniden üretildi.

## Kaynaklar

- Grantham uzaklığı: [Grantham, Science 1974](https://doi.org/10.1126/science.185.4154.862)
- BLOSUM62: [Henikoff ve Henikoff, PNAS 1992](https://doi.org/10.1073/pnas.89.22.10915)
- Hidropati ölçeği: [Kyte ve Doolittle, J Mol Biol 1982](https://doi.org/10.1016/0022-2836(82)90515-0)

Sabit tablolar çalışma zamanında internetten çekilmez; sürüm kontrollü özellik
üretici içinde tutulur ve boyut/simetri testleriyle denetlenir.

## Deney düzeni

İki lojistik regresyon aynı 257.035 gene-disjoint train ve 55.068 validation
kaydı, `seed=2027` ve sabit `0.5` eşiğiyle karşılaştırıldı:

1. Tam HGVS referans profili.
2. Aynı profil + missense biyolojik özellikleri.

55.093 satırlık test bölümü geliştirme veri çerçevesine yüklenmedi. Model veya
eşik için tarama yapılmadı.

## Genel validation sonucu

| Profil | ROC-AUC | AP | Precision | Recall | Özgüllük | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| HGVS referans | 0.9824 | 0.9810 | 0.9505 | 0.8902 | 0.9568 | 0.9194 |
| + missense biyolojik özellikler | 0.9859 | 0.9845 | 0.9558 | 0.9210 | 0.9603 | 0.9381 |

Zenginleştirilmiş modelde yanlış negatif sayısı 2.919'dan 2.100'e, yanlış
pozitif sayısı 1.231'den 1.132'ye düşmüştür.

## Missense validation sonucu

| Profil | n | Benign | Pathogenic | ROC-AUC | AP | Precision | Recall | Özgüllük | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| HGVS referans | 7.075 | 3.884 | 3.191 | 0.5696 | 0.5078 | 0.5262 | 0.3149 | 0.7670 | 0.3940 |
| + missense biyolojik özellikler | 7.075 | 3.884 | 3.191 | 0.7513 | 0.6878 | 0.6931 | 0.5704 | 0.7925 | 0.6258 |

Zenginleştirilmiş eksi referans farkları:

- ROC-AUC: `+0.181659`
- AP: `+0.179989`
- precision: `+0.166891`
- recall: `+0.255406`
- özgüllük: `+0.025489`
- F1: `+0.231712`

Missense yanlış negatifleri 2.186'dan 1.371'e, yanlış pozitifleri 905'ten
806'ya düşmüştür.

## Yorum ve sınırlar

Kazancın genel sonuçtan çok missense alt grubunda yoğunlaşması, özelliklerin
amaçlanan biyolojik boşluğu kısmen doldurduğunu destekler. Ancak bu sonuç tek
sabit validation ölçümüdür. Yönlü aminoasit çiftleri ClinVar'daki gönderim ve
etiketleme örüntülerini de öğrenebilir; katsayılar nedensel biyolojik etki olarak
yorumlanmamalıdır.

Protein konumu, korunmuşluk, yapısal bölge ve dış popülasyon frekansı henüz
kullanılmamıştır. Yeni profili ana model adayı yapmadan önce train-only
gene-disjoint çapraz doğrulamada fold kararlılığı ölçülmeli, ardından geliştirme
sürecinde kullanılmamış dış veriyle doğrulama planlanmalıdır.

## Yeniden çalıştırma

```powershell
.\.venv\Scripts\python.exe scripts\build_missense_features.py
.\.venv\Scripts\python.exe scripts\run_missense_feature_ablation.py
```
