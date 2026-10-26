# ClinGen kaydedilmiş tahminlerinde hata analizi

Tarih: 26 Ekim 2026

## Sonuç

Birincil değerlendirmedeki 80 yanlış negatifin 72'si missense grubundadır.
Hatalar, eğitimde hiç görülmeyen aminoasit değişimleriyle açıklanamıyor: 187
missense kaydının tamamındaki yönlü aminoasit çifti modelin fit verisinde
görülmüştür. Kategorik model alanlarında bilinmeyen değer sayısı da 0'dır.

Bu, dağılımın tamamen aynı olduğunu veya tek bir hata nedenini kanıtlamaz.
Veri seçimi, sınıf dağılımı, anotasyon kaynağı ve hastalık bağlamı farklıdır.
Analiz sonuç sonrası açıklayıcı incelemedir; yeni bir doğrulama deneyi değildir.

## Sabit kalanlar ve kullanılan girdiler

- Yalnızca kaydedilmiş olasılık/kararlar okundu; model yüklenmedi veya eğitilmedi.
- Birincil değerlendirme ve tahmin dosyalarının SHA-256 değerleri öncesinde ve
  sonrasında doğrulandı; sabit `0.50` eşik korundu.
- Eğitim desteği yalnızca 312.103 train+validation kaydının özelliklerinden
  sayıldı. İç test hedefleri okunmadı.
- Gen, kaynak, protein etkisi ve aminoasit çifti gruplarında FN/Pathogenic ve
  FP/Benign oranları hesaplandı. Paydası olmayan oranlar boş bırakıldı.
- CSV'de Wilson %95 aralıkları ve sınıf başına asgari 5 kayıt bayrakları var.
  Bu aralıklar satır bağımsızlığı varsayar; gen kümelenmesi ve çoklu inceleme
  için düzeltilmemiştir. Asgari 5 kayıt istatistiksel yeterlilik garantisi değildir.

## Hata yükünün yoğunlaştığı gruplar

| Protein etkisi | Pathogenic | Benign | FN | FP | FN/Pathogenic |
| --- | ---: | ---: | ---: | ---: | ---: |
| Missense | 148 | 39 | 72 | 10 | %48,65 |
| Protein anotasyonu yok | 15 | 17 | 6 | 0 | %40,00 |
| Synonymous | 2 | 8 | 2 | 0 | %100,00 |
| Stop-gained/changed | 22 | 1 | 0 | 1 | %0,00 |
| Frameshift | 9 | 0 | 0 | 0 | %0,00 |

Synonymous satırında yalnızca iki Pathogenic örnek vardır; bu oran genel bir
synonymous performans tahmini olarak kullanılmamalıdır.

### Genler

| Gen | Pathogenic | Benign | FN | FP | FN/Pathogenic |
| --- | ---: | ---: | ---: | ---: | ---: |
| RYR1 | 23 | 0 | 16 | 0 | %69,57 |
| TP53 | 34 | 1 | 10 | 0 | %29,41 |
| ABCA4 | 14 | 1 | 5 | 0 | %35,71 |
| PTEN | 11 | 0 | 4 | 0 | %36,36 |
| OTC | 6 | 9 | 2 | 2 | %33,33 |

RYR1 ve TP53 birlikte 26/80 yanlış negatif içerir (`%32,5`). Bunlar hata yükü
sayımlarıdır. Özellikle Benign örneği olmayan RYR1 veya PTEN için gen bazlı
özgüllük/ROC karşılaştırması yapılamaz. Genler burada inceleme alanıdır; model
girdisine gen kimliği eklenmedi.

### Missense aminoasit çiftleri

| Çift | Pathogenic | Benign | FN | FP | FN/Pathogenic |
| --- | ---: | ---: | ---: | ---: | ---: |
| R>C | 10 | 1 | 8 | 1 | %80,00 |
| R>H | 8 | 2 | 8 | 0 | %100,00 |
| E>K | 8 | 1 | 7 | 1 | %87,50 |
| R>Q | 6 | 0 | 6 | 0 | %100,00 |
| V>M | 5 | 2 | 5 | 0 | %100,00 |
| D>N | 5 | 2 | 5 | 0 | %100,00 |

Bu altı çift, 72 missense yanlış negatifin 39'unu (`%54,17`) içerir. Küçük
örnek sayıları ve sonuç sonrası seçim nedeniyle bunlar keşifsel inceleme
öncelikleridir; bu çiftlere özel yeni karar kuralları oluşturulmadı.

176 missense kaydının aminoasit çifti fit verisinde en az 100 kez görülmüştür;
bu kayıtlarda 70 yanlış negatif vardır. Diğer 11 kaydın çifti 20–99 kez
görülmüştür. Dolayısıyla hatalar yalnızca nadir veya bilinmeyen çiftlerde
toplanmıyor; kategori görülmüşlüğü bağlamsal genellemeyi garanti etmiyor.

## Anotasyon kaynağı

| Kaynak / kapsam | Pathogenic | Benign | FN | FP | FN/Pathogenic |
| --- | ---: | ---: | ---: | ---: | ---: |
| Eski ClinVar Name / tümü | 96 | 58 | 45 | 11 | %46,88 |
| eRepo başlığı / tümü | 101 | 8 | 35 | 1 | %34,65 |
| Eski ClinVar Name / missense | 76 | 33 | 39 | 10 | %51,32 |
| eRepo başlığı / missense | 72 | 6 | 33 | 0 | %45,83 |

Kaynaklar aynı örneklerin alternatif anotasyonları değildir; gen, sonuç türü ve
etiket dağılımları farklıdır. Bu tablo, kaynak değiştirmenin performansı
artıracağını kanıtlamaz. Tek Registry başlığı kaydı bu karşılaştırmaya yetmez.

## Özellik gösteriminin sınırı

5 missense grubu, **modelin 25 girdi alanında birebir aynı değerleri** taşıdığı
hâlde farklı uzman etiketleri içerir. Bu gruplarda 12 varyant, 6 FN ve 1 FP
vardır; kaydedilmiş olasılıklar grup içinde aynıdır.

Örneğin R>H çifti taşıyan OTC ve SLC6A8 kayıtları aynı girdilerle `0.4471592689`
skoru almış, fakat biri Pathogenic diğeri Benign olarak etiketlenmiştir. Mevcut
gösterimi kullanan deterministik bir model bu iki kaydı ayıramaz. Bu gözlem
bağlam bilgisinin eksikliğini somutlaştırır; 80 FN'nin tamamını açıklamaz ve
hangi ek özelliğin işe yarayacağını henüz göstermez.

## Uzman kanıtından üç inceleme örneği

Aşağıdaki açıklamalar aynı tarihli eRepo snapshot'ındaki `summaryDesc`
alanlarının özetidir; yeni klinik yorum yapılmadı. Tam metinler ve yorum UUID'leri
yerel `error_curation_evidence.tsv` dosyasında saklanır.

| VariationID | Kaydedilmiş skor | Hata | Snapshot uzman özetindeki ilgili bağlam |
| --- | ---: | --- | --- |
| 569548 / ACADVL | 0.0107403065 | FN | Synonymous gösterime rağmen RNA/ekzon atlama kanıtı belirtilir. |
| 38682 / HBB | 0.0075507278 | FN | Synonymous gösterime rağmen splice öngörüsü ve klinik/segregasyon kanıtları birlikte tartışılır. |
| 190351 / GFAP | 0.9975222204 | FP | Stop-gained gösterimin tek başına yeterli olmadığı; Alexander hastalığı bağlamında mekanizma ve frekans kanıtıyla Benign yorumu verildiği açıklanır. |

Mevcut kaba `protein_effect` ve HGVS bölge kategorileri, bu kanıtları ayrıntılı
biçimde kodlamaz. Bunlar uzman denetimi ve sonraki özellik araştırması için
örneklerdir; bu görülmüş kohort üzerinde model veya eşik ayarlanmadı.

## Çıktılar ve sonraki çalışma

```powershell
.\.venv\Scripts\python.exe scripts\analyze_clingen_errors.py
```

Girdi kilidi: `config/clingen_error_analysis_inputs.json`.
Yerel çıktılar `reports/generated/clingen_error_analysis/` altındadır:

- `subgroup_error_rates.csv`: tüm grup sayıları, oranları ve Wilson aralıkları.
- `saved_error_examples.tsv`: 92 hatanın kimlik, Name, skor ve özellik bağlamı.
- `missense_feature_collisions.tsv`: aynı girdilerle farklı etiket taşıyan gruplar.
- `error_curation_evidence.tsv`: hata kayıtlarının kaynak uzman özetleri.

Sonraki çalışma, mevcut modele dokunmadan splice ve fonksiyonel bağlam
anotasyonlarının kapsanabilirliğini araştırmak ve yeni deney için kullanılmamış
bir değerlendirme kohortu ayırmaktır. Ek anotasyonlar label/kanıt özetlerinden
hedef sızıntısı yaratacak şekilde türetilmemelidir. Başarı ancak ayrı geliştirme
deneyi ve yeni doğrulama verisiyle gösterilebilir.
