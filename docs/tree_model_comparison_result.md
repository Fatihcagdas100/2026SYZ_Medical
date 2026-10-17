# Ağaç modeli karşılaştırma sonucu

Tarih: 18 Ekim 2026

## Deney düzeni

Lojistik regresyon ve `HistGradientBoostingClassifier` şu koşullarda
karşılaştırıldı:

- aynı 257.035 gene-disjoint train kaydı,
- aynı 55.068 gene-disjoint validation kaydı,
- aynı tam HGVS özellik seti,
- açık VCF eksiklik bayrağı olmadan,
- `seed=2027`,
- sabit `0.5` karar eşiği.

55.093 satırlık test bölümü geliştirme veri çerçevesine yüklenmedi ve test
metriği üretilmedi. Ağaç modeli için tek bir yapılandırma koşudan önce
sabitlendi; validation üzerinde hiperparametre taraması yapılmadı:

- `learning_rate=0.08`
- `max_iter=200`
- `max_leaf_nodes=31`
- `min_samples_leaf=50`
- `l2_regularization=1.0`
- `early_stopping=False`

Bu deney bir model seçimi taraması değil, tek ve yeniden üretilebilir bir ağaç
baseline karşılaştırmasıdır.

## Genel validation sonucu

| Model | ROC-AUC | AP | Precision | Recall | Özgüllük | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Lojistik regresyon | 0.9824 | 0.9810 | 0.9502 | 0.8908 | 0.9564 | 0.9195 |
| HistGradientBoosting | 0.9839 | 0.9821 | 0.9506 | 0.8984 | 0.9564 | 0.9238 |

Ağaç eksi lojistik farkları:

- ROC-AUC: `+0.001519`
- AP: `+0.001124`
- recall: `+0.007602`
- F1: `+0.004223`
- özgüllük: `0.000000`

Yanlış negatif sayısı 2.902'den 2.700'e düşmüştür. Yanlış pozitif sayısı iki
modelde de 1.241'dir.

## Kritik alt gruplar

| Grup | Model | ROC-AUC | AP | Precision | Recall | Özgüllük |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Missense | Lojistik | 0.5703 | 0.5089 | 0.5263 | 0.3196 | 0.7636 |
| Missense | Ağaç | 0.5574 | 0.4998 | 0.5159 | 0.3149 | 0.7572 |
| Synonymous | Lojistik | 0.5422 | 0.0118 | 0.0000 | 0.0000 | 1.0000 |
| Synonymous | Ağaç | 0.5775 | 0.0115 | 0.0000 | 0.0000 | 1.0000 |
| Frameshift | Lojistik | 0.5615 | 0.9949 | 0.9940 | 0.9997 | 0.0000 |
| Frameshift | Ağaç | 0.5137 | 0.9945 | 0.9940 | 1.0000 | 0.0000 |
| Stop-gained/changed | Lojistik | 0.4667 | 0.9967 | 0.9969 | 1.0000 | 0.0000 |
| Stop-gained/changed | Ağaç | 0.3782 | 0.9959 | 0.9969 | 1.0000 | 0.0000 |

Missense grubunda ağaç modeli lojistiğe göre daha kötüdür:

- ROC-AUC farkı: `-0.012884`
- AP farkı: `-0.009084`
- recall farkı: `-0.004701`
- F1 farkı: `-0.006612`

Synonymous gruptaki 46 patojenik kaydın hiçbirini iki model de `0.5` eşiğinde
yakalamamıştır. Frameshift ve stop-gained/changed gruplarında iki model de bütün
benign kayıtları patojenik tahmin ettiği için özgüllük sıfırdır.

## Anotasyon varlığına göre sonuç

| Grup | Model | ROC-AUC | AP | Recall | Özgüllük |
| --- | --- | ---: | ---: | ---: | ---: |
| Protein anotasyonu var | Lojistik | 0.9659 | 0.9844 | 0.8999 | 0.8829 |
| Protein anotasyonu var | Ağaç | 0.9653 | 0.9839 | 0.9014 | 0.8794 |
| Protein anotasyonu yok | Lojistik | 0.9750 | 0.9405 | 0.8367 | 0.9929 |
| Protein anotasyonu yok | Ağaç | 0.9804 | 0.9584 | 0.8807 | 0.9946 |

Genel artışın önemli bölümü protein anotasyonu olmayan gruptan gelmektedir. Bu
grupta ağaç modeli yanlış negatifleri 623'ten 455'e, yanlış pozitifleri 135'ten
102'ye düşürmüştür. Protein anotasyonu bulunan grupta ROC-AUC ve AP küçük
ölçüde azalmıştır.

## Karar

Ağaç modeli genel validation metriklerinde küçük bir artış sağlamıştır; fakat
asıl zor grup olan missense varyantlarını iyileştirmemiştir. Dolayısıyla yalnızca
genel ROC-AUC'ye bakarak ağaç modelini yeni ana model ilan etmek doğru değildir.

Sonuç, hataların model ailesinden çok mevcut özelliklerin bilgi sınırında
biriktiği hipoteziyle uyumludur. Sıradaki geliştirme adımı daha büyük veya daha
çok ayarlanmış bir model değil; missense varyantlar içinde ayrım sağlayabilecek
hedefli biyolojik özelliklerin sızıntı ve lisans koşullarıyla birlikte
değerlendirilmesidir.

## Yeniden çalıştırma

```powershell
.\.venv\Scripts\python.exe scripts\compare_tree_model.py
```

Yerel ayrıntılı çıktılar:

- `reports/generated/tree_model_comparison/validation_model_comparison.csv`
- `reports/generated/tree_model_comparison/validation_subgroup_comparison.csv`
- `reports/generated/tree_model_comparison/validation_metric_deltas.csv`
- `reports/generated/tree_model_comparison/tree_model_comparison.md`
