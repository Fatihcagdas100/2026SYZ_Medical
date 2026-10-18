# Missense özellikleri: train-only gene-disjoint OOF kararlılık kontrolü

Tarih: 18 Ekim 2026

## Amaç ve veri sınırı

Bu deney, missense biyolojik özelliklerinin sabit validation'daki kazanımının tek
bir gene-disjoint ayrımına özgü olup olmadığını inceler. Yalnızca 257.035 satırlı
`train` bölümü kullanıldı. Sabit 55.068 satırlı validation bölümü ve 55.093
satırlı kilitli test bölümü yüklenmedi.

İki lojistik profil, aynı `seed=2027`, sabit `0.5` eşik ve 5 katlı
`StratifiedGroupKFold` ile çalıştırıldı:

1. Tam HGVS referans profili.
2. Aynı profil + aminoasit ikame biyolojisi özellikleri.

Her gen bileşeni yalnızca bir holdout fold'unda yer aldı. İki profilin OOF
tahminlerindeki `VariationID` kümesi ve fold ataması birebir aynı olarak
denetlendi.

## Birleştirilmiş OOF sonuçları

| Grup | Profil | n | ROC-AUC | AP | Precision | Recall | Özgüllük | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Tüm varyantlar | Tam HGVS referans | 257.035 | 0.9829 | 0.9823 | 0.9471 | 0.8959 | 0.9533 | 0.9208 |
| Tüm varyantlar | HGVS + missense biyoloji | 257.035 | 0.9862 | 0.9857 | 0.9558 | 0.9262 | 0.9600 | 0.9407 |
| Missense | Tam HGVS referans | 33.744 | 0.5503 | 0.4939 | 0.4958 | 0.3200 | 0.7372 | 0.3889 |
| Missense | HGVS + missense biyoloji | 33.744 | 0.7358 | 0.6811 | 0.6795 | 0.5691 | 0.7832 | 0.6194 |

Missense grubunda zenginleştirilmiş profilin farkları şöyledir:

- ROC-AUC: `+0.185489`
- AP: `+0.187118`
- precision: `+0.183710`
- recall: `+0.249138`
- özgüllük: `+0.046015`
- F1: `+0.230489`

Fold bazında genel ROC-AUC ortalaması referans için `0.9829 ± 0.0040`,
zenginleştirilmiş profil için `0.9862 ± 0.0029` olmuştur.

## Yorum

Missense kazancının hem sabit validation'da hem de yalnızca eğitim verisindeki
gen-korumalı OOF tahminlerinde aynı yön ve benzer büyüklükte görülmesi, bulgunun
tek bir geliştirme ayrımının rastlantısı olmadığına dair daha güçlü bir iç kanıt
sağlar. Bu deney yine ClinVar tabanlı geliştirme verisindedir: yönlü aminoasit
çiftleri veri kaynağının seçilme veya etiketleme örüntülerini de taşıyabilir.
Sonuç dış veri genellenmesi, klinik tanı performansı veya nedensel mekanizma
kanıtı olarak yorumlanmamalıdır.

## Yeniden çalıştırma

```powershell
.\.venv\Scripts\python.exe scripts\run_missense_gene_disjoint_cv.py
```

Üretilen yerel raporlar `reports/generated/missense_gene_disjoint_cv/` altında
bulunur. Bu klasör ham deney çıktısıdır ve sürüm kontrolüne alınmaz.
