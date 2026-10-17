# İç gene-disjoint değerlendirme sonucu

Tarih: 17 Ekim 2026

## Sonuç

Önceden dondurulmuş model, özellik seti ve `0.3623531848` karar eşiğiyle
55.093 varyant ve son modelin eğitiminde görülmemiş 1.798 gen grubundan oluşan
gene-disjoint test bölümü tek sefer değerlendirilmiştir. Önceden belirlenen beş
proje-içi başarı kapısının tamamı geçilmiştir.

Ayrım geçiş denetimi, bu testteki 46.830 varyantın (`%85,00`) daha önceki
varyant-düzeyi deneylerde train veya validation rolünde bulunduğunu daha sonra
göstermiştir. Dolayısıyla doğrudan son-model train/test sızıntısı olmasa da bu
sonuç geliştirme süreci açısından tamamen görülmemiş bağımsız nihai test olarak
sunulamaz. Doğru adlandırma **iç gene-disjoint değerlendirme**dir.

| Test metriği | Değer |
| --- | ---: |
| ROC-AUC | 0.9872 |
| Average precision | 0.9870 |
| Log loss | 0.1356 |
| Brier skoru | 0.0406 |
| Duyarlılık | 0.9661 |
| Özgüllük | 0.8762 |
| Precision | 0.8793 |
| F1 | 0.9206 |
| Dengeli doğruluk | 0.9212 |
| Matthews korelasyon katsayısı | 0.8434 |

Donmuş eşikte karışıklık matrisi:

|  | Tahmin: Benign | Tahmin: Pathogenic |
| --- | ---: | ---: |
| Gerçek: Benign | 24.974 | 3.527 |
| Gerçek: Pathogenic | 902 | 25.690 |

## Önceden belirlenmiş proje-içi başarı kapıları

| Başarı kapısı | Gözlenen | Sınır | Sonuç |
| --- | ---: | ---: | --- |
| ROC-AUC | 0.9872 | ≥ 0.97 | Geçti |
| Average precision | 0.9870 | ≥ 0.97 | Geçti |
| Duyarlılık | 0.9661 | ≥ 0.93 | Geçti |
| Özgüllük | 0.8762 | ≥ 0.85 | Geçti |
| Validation'a göre ROC-AUC düşüşü | -0.0048 | ≤ 0.02 | Geçti |

Negatif ROC-AUC düşüşü, test ROC-AUC değerinin validation değerinden yaklaşık
0.0048 daha yüksek olduğunu gösterir.

Bu kapılar test açılmadan önce belirlenmiştir; ancak klinik geçerlilik veya
TEKNOFEST tarafından tanımlanmış başarı standartları değildir.

## Teknik dtype düzeltme denetimi

İlk koşu sırasında `chromosome` sütunu için pandas karışık veri tipi uyarısı
üretmiştir. Eğitim hattında kategorik alanların string olarak okunduğu
doğrulandığı için orijinal koşu korunmuş, model veya eşik değiştirilmeden ayrı
bir teknik düzeltme denetimi yapılmıştır.

- String veri tipi düzeltmesinden etkilenen olasılık: 208
- Sınıf kararı değişen kayıt: 0
- ROC-AUC farkı: -0.000005
- Average precision farkı: -0.000002
- En büyük tekil olasılık farkı: 0.033284
- Düzeltme sonrasında başarı kapıları: tamamı geçti

Bu belgede yukarıdaki ana tabloda dtype-düzeltilmiş değerler raporlanmıştır.
Orijinal koşu ve tahmin dosyaları silinmemiş veya üzerlerine yazılmamıştır.
Sınıf kararlarının değişmemesi ROC-AUC/AP gibi eşikten bağımsız ölçülerin tam
olarak aynı olduğu anlamına gelmez: ROC-AUC `0.987192`den `0.987187`ye, AP
`0.986995`ten `0.986992`ye değişmiştir.

## Ayrım geçiş denetimi

| Önceki varyant-düzeyi rol | Son gene-disjoint testteki varyant | Oran |
| --- | ---: | ---: |
| Train | 38.630 | %70,12 |
| Validation | 8.200 | %14,88 |
| Test | 8.263 | %15,00 |

Önceki train ve validation toplamı 46.830 varyanttır (`%85,00`). Son modelin
parametreleri bu eski ayrımla eğitilmemiştir; ancak önceki deney sonuçları hata
analizini ve HGVS özellik geliştirmesini yönlendirmiştir. Bu nedenle testin
süreç düzeyindeki bağımsızlığı sınırlıdır.

## Bütünlük bilgileri

- Protokol kilidi SHA-256:
  `7497977afc825922e877fc24c333d47502f0bc7fc5e60f0751eaaffdfdf901f9`
- Model artifact SHA-256:
  `a4de323fbc301f9ff0fa28e30148648cce547343cb3c68491c9d394470d00c4c`
- Özellik tablosu SHA-256:
  `4c0814b5b0bf4d47ff17357aa8870d866ab61846f2ce1ccd625a078e1a5f3e3a`
- Gene-disjoint ayrım SHA-256:
  `7c115e25dd1c7edc1bc566ad12c14e829aa02252741dda8d080f4ac3cc796aa7`

## Yorum ve sınırlar

Bu sonuç modelin aynı ClinVar veri kaynağı içinde son model eğitiminde
görülmemiş gen gruplarında güçlü performans gösterdiğini destekler. Bununla
birlikte önceki geliştirme maruziyeti nedeniyle süreçten bağımsız bir test veya
bağımsız dış veri doğrulaması değildir. Model klinik tanı sistemi olarak
kullanılamaz; araştırma amaçlı varyant önceliklendirme prototipi olarak ve uzman
değerlendirmesiyle birlikte ele alınmalıdır.
