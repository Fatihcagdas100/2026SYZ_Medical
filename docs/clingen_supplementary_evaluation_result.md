# ClinGen uzman-kürasyonlu ek değerlendirme sonucu

Tarih: 26 Ekim 2026

Kaydedilmiş tahminlerin takip hata analizi:
`docs/clingen_error_analysis_result.md`. Birincil ölçüm ve aşağıdaki sonuçlar
takip analizinde değiştirilmedi.

## Birincil sonuç

Dondurulmuş `hgvs_missense_development_refit` modeli, 264 varyantlık kabul edilmiş
ClinGen kohortunda bir kez ölçüldü. Karar eşiği önceden sabitlenmiş `0.50` olarak
korundu. Genel ROC-AUC **0.7932**, missense ROC-AUC **0.6642** bulundu. Genel
duyarlılık `%59,60`; 198 Pathogenic kaydın 80'i yanlış negatifti.

Bu sonuç iç geliştirme performansının bu seçilmiş kohorta aynı düzeyde
taşınmadığını gösterir. Ölçüm sonrası model, özellikler veya eşik değiştirilmedi.
Önceden bu kohorta özgü sayısal başarı kapıları tanımlanmadığı için sonuca
sonradan 'geçti/geçmedi' eşiği atanmadı.

## Kohort ve ölçüm düzeni

- Kapsam: `expert_curated_supplementary`; kaynak bağımsızlığı doğrulanmış değil.
- Tüm kohort: 264 varyant, 198 Pathogenic ve 66 Benign.
- Missense: 187 varyant, 148 Pathogenic ve 39 Benign.
- Yeni/görülmemiş gen kaydı: 0; bu alt grup için performans hesaplanmadı.
- Protein anotasyonu olmayan kayıt: 32.
- Bootstrap: 1.000 tekrar, `seed=2027`, sınıf sayıları korunan satır örneklemesi.
- Ölçüm ayarları: `config/clingen_supplementary_evaluation.json`.

Güven aralıkları sınıf prevalansını sabit tutar; aynı gen içindeki kayıtların
bağımlılığını modellemez. Küçük Benign alt sınıfları nedeniyle belirsizlik
aralıkları geniştir. Bu aralıklar yöntem kapsamında yorumlanmalıdır.

## Metrikler

| Metrik | Tüm kohort | %95 GA | Missense | %95 GA |
| --- | ---: | --- | ---: | --- |
| ROC-AUC | 0.7932 | 0.7357–0.8548 | 0.6642 | 0.5773–0.7546 |
| Average precision | 0.9111 | 0.8761–0.9436 | 0.8754 | 0.8289–0.9213 |
| Precision | 0.9077 | 0.8640–0.9520 | 0.8837 | 0.8313–0.9383 |
| Duyarlılık | 0.5960 | 0.5303–0.6667 | 0.5135 | 0.4390–0.5878 |
| Özgüllük | 0.8182 | 0.7273–0.9091 | 0.7436 | 0.6154–0.8718 |
| F1 | 0.7195 | 0.6646–0.7733 | 0.6496 | 0.5804–0.7137 |

Pathogenic oranı genel kohortta `%75`, missense grubunda `%79,14` olduğundan
yüksek AP tek başına güçlü ayırt etme veya yeterli duyarlılık anlamına gelmez.
Prevalans düzeyleri, AP için ilgili taban karşılaştırmalarıdır.

## Karar hataları

| Grup | TN | FP | FN | TP |
| --- | ---: | ---: | ---: | ---: |
| Tüm kohort | 54 | 12 | 80 | 118 |
| Missense | 29 | 10 | 72 | 76 |

Kaydedilmiş tahminler üzerinden hata sayımları ayrıca doğrulandı. 80 yanlış
negatifin 72'si (`%90`) missense, 6'sı protein anotasyonsuz ve 2'si synonymous
olarak ayrışır. 12 yanlış pozitifin 10'u missense, 1'i stop-gained/changed ve
1'i start-lost/uncertain grubundadır. Bu, hata yükünün missense grubunda
yoğunlaştığını gösterir; tek başına hatanın biyolojik veya teknik nedenini
kanıtlamaz.

## İç sonuçlarla ilişki ve yorum sınırı

Önceki train-only OOF'da genel ROC-AUC `0.9862`, missense ROC-AUC `0.7358`
bulunmuştu. Buradaki kohort farklı seçilme, etiket ve anotasyon kaynakları
taşır; varyantların çoğu missense'dir. Bu iki deneyin farkını aynı dağılımda
kontrollü bir düşüş veya tek bir nedene bağlanan etki olarak yorumlamıyoruz.

Mevcut bulgu, özellikle missense ve sabit eşikte duyarlılık açısından
genellenebilirliğin sınırlı kaldığına işaret eder. Kaynak ilişkileri ve tüm
genlerin önceden görülmüş olması nedeniyle bu sonuç bağımsız klinik dış
doğrulama olarak sunulamaz. Sonraki geliştirmeler bu kohort artık görülmüş
olduğu için başka, kullanılmamış bir doğrulama kohortuyla teyit edilmelidir.

## Bütünlük ve tekrar kuralı

Etiketler ölçüm için okunmadan önce `EVALUATION_STARTED.json` ve aday model
klasöründe ham kohort hash'ine bağlı registry kaydı oluşturuldu. Başarılı
ölçüm manifesti ve kayıtlı tahminler saklandı; yeniden değerlendirme yapılmadı.
Eski iç testte yeni ölçüm yapılmadı.

- Model artifact SHA-256:
  `315b8a19986dcb0e20c010066c24f0baae62255d042828aefcc1fd78905139f7`
- Özellik dosyası SHA-256:
  `a971ad3819e46bec396735211230e398f9dce3c1d1e324b7a0d014b2232a4683`
- Sonuç manifesti SHA-256:
  `d4c72dce69278d1343d74158ca8e8dfe08fc2b71f11e7b02ce824c60825feccd`
- Kaydedilmiş tahminler SHA-256:
  `291c244b4a198e400ab87030a0e4bbb7381ea2cfd15c0a2ef5c1da619da6bafb`
- Ölçüm ayar dosyası SHA-256:
  `96a8351eee23f053c3f4b586dabc6b99bc29f0f0c7a6037364e01662b50eb18d`

Yerel ham sonuçlar `reports/generated/clingen_supplementary_evaluation/`
altındadır. Bu belge sürüm kontrolündeki sonuç özetidir; komutu tekrar
çalıştırmak yerine kayıtlı manifest ve tahminler incelenmelidir.
