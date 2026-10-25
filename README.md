# 2026SYZ

TEKNOFEST Sağlıkta Yapay Zeka 2027 hazırlığı kapsamında genetik varyantların
patojenik/benign sınıflandırılması üzerine yürütülen çalışma.

## Mevcut durum

- Katı `GRCh38 + germline + tam Pathogenic/Benign` politikasıyla 367.196
  benzersiz ClinVar varyantı hazırlandı.
- Kesin genom koordinatları ve gen kimliği modelden dışlandı; kromozom yalnızca
  kategorik özellik olarak kullanıldı.
- Son lojistik model HGVS sonuç özellikleriyle gene-disjoint train bölümünde
  eğitildi. Train-only 5-fold gene-disjoint ROC-AUC `0.9828 ± 0.0041` oldu.
- Gene-disjoint iç test 55.093 varyant ve 1.798 model-eğitiminde görülmemiş gen
  grubunda bir kez değerlendirildi; ROC-AUC `0.9872`, AP `0.9870`, duyarlılık
  `%96,61` ve özgüllük `%87,62` bulundu.
- Ayrım geçiş denetimi, bu iç test varyantlarının `%85,00`'ının daha önceki
  varyant-düzeyi deneylerde train veya validation rolünde bulunduğunu gösterdi.
  Bu nedenle sonuç tamamen süreçten bağımsız nihai test değil, **iç
  gene-disjoint değerlendirme** olarak raporlanır.
- Sonraki doğrulama adımı geliştirme sürecinde hiç kullanılmamış bağımsız bir
  dış veri kaynağıdır. `HGVS + missense` dış doğrulama adayı yalnızca train ve
  validation ile donduruldu; dış kohort henüz değerlendirilmedi. Model klinik
  kullanıma hazır değildir.

## Python ortamı

VS Code için yorumlayıcı:

```text
C:\Users\Asus\Desktop\2026SYZ\.venv\Scripts\python.exe
```

Terminalden doğrulama:

```powershell
.\.venv\Scripts\python.exe --version
```

Modelleme bağımlılıklarını proje ortamına kurmak için:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Ortamı doğrulanan sürümlerle birebir yeniden kurmak için:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
```

`requirements.txt` doğrudan kullandığımız paketleri, `requirements-lock.txt`
ise yeniden üretim için kurulan bütün doğrudan ve dolaylı paket sürümlerini
sabitler.

## Veri

Ham veri büyük olduğu ve harici bir kaynaktan yeniden indirilebildiği için Git
repository'sine eklenmez. Beklenen yerel konum:

```text
data/raw/variant_summary.txt.gz
```

Kaynak: [NCBI ClinVar indirme sayfası](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/)

Doğrudan dosya: [variant_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)

## Veri denetimi

Tam dosyayı belleğe yüklemeden denetlemek için:

```powershell
.\.venv\Scripts\python.exe scripts\audit_variant_summary.py
```

Önce hızlı bir deneme yapmak için:

```powershell
.\.venv\Scripts\python.exe scripts\audit_variant_summary.py --max-rows 10000
```

Betik terminalde bütün sütun adlarını ve seçili alanlardan örnek kayıtları
gösterir. Ayrıntılı çıktılar `reports/generated/clinvar_audit/` altında üretilir:

- `summary.json`: tarama ve tekrar özeti,
- `columns.csv`: sıralı sütun listesi,
- `column_profile.csv`: boş, `-`, `na` ve dolu değer sayıları,
- `categorical_distributions.csv`: seçili sütunların değer dağılımları,
- `sample_rows.tsv`: ilk örnek kayıtlar ve 43 sütunun tamamı,
- `variation_id_repeat_examples.csv`: en sık tekrarlanan kimlik örnekleri.

`VariationID` tekrarları geçici bir SQLite dosyasında izlendiği için milyonlarca
kimlik bellekte tutulmaz. Geçici indeks işlem sonunda otomatik olarak silinir.

### İlk ikili alt küme profilini çıkarma

Henüz yeni veri seti yazmadan `GRCh38 + germline + tam Pathogenic/Benign`
filtrelerinin kaç kayıt bıraktığını görmek için:

```powershell
.\.venv\Scripts\python.exe scripts\profile_initial_binary_subset.py
```

Sonuçlar `reports/generated/initial_binary_subset/` altında oluşturulur. Bu
adımda `Likely pathogenic`, `Likely benign`, birleşik, belirsiz ve çelişkili
etiketler bilinçli olarak ikili hedefe katılmaz.

### İlk türetilmiş veri setini oluşturma

Filtre profilindeki katı politikayla her `VariationID` için tek satır üretmek
için:

```powershell
.\.venv\Scripts\python.exe scripts\build_initial_binary_dataset.py
```

Çıktı `data/processed/clinvar_grch38_germline_exact_binary.tsv.gz` konumuna
yazılır. `Pathogenic=1`, `Benign=0` eşlemesi kullanılır. Birden fazla farklı
kromozom/koordinat temsili bulunan `VariationID` değerleri ilk sade veri
setinden dışlanır. Kaynak ve çıktı SHA-256 değerleri, filtre politikası ve sütun
kullanım uyarıları `reports/generated/initial_binary_dataset/manifest.json`
içinde kaydedilir.

### Temel özellik tablosunu oluşturma

Kimlik, gen ve kesin genom koordinatları gibi ezberleme riski yüksek alanları
dışarıda tutarak ilk yapısal özellik tablosunu üretmek için:

```powershell
.\.venv\Scripts\python.exe scripts\build_baseline_features.py
```

Çıktı `data/processed/clinvar_baseline_features.tsv.gz` konumuna yazılır. Bu
dosyada varyant türü, kromozom, allel uzunlukları, uzunluk farkı, kanonik SNV
durumu, baz değişimi ve mutasyon sınıfı bulunur. `VariationID` yalnızca satır
izlenebilirliği içindir; model girdisi değildir.

### Sabit train/validation/test ayrımı

Sınıf oranlarını koruyan deterministik `%70 / %15 / %15` ayrımı oluşturmak
için:

```powershell
.\.venv\Scripts\python.exe scripts\build_stratified_split.py
```

Ayrım `seed=2027` ve SHA-256 sıralamasıyla üretilir. Atamalar
`data/processed/baseline_stratified_split.tsv`, yöntem ve sayımlar ise
`reports/generated/baseline_split/manifest.json` altında tutulur. Bu ilk ayrım
varyant düzeyindedir; görülmemiş genlere genelleme deneyi değildir.

### İlk lojistik regresyon baseline'ı

Modeli yalnızca `train` bölümünde eğitip sadece `validation` bölümünde
değerlendirmek için:

```powershell
.\.venv\Scripts\python.exe scripts\train_logistic_baseline.py
```

Betik sabit `0.5` karar eşiğiyle iki kontrollü model kurar: ana model açık VCF
eksiklik göstergesini kullanmaz; karşılaştırma modeli bu göstergeyi ekler.
`VariationID` model girdisi değildir. `test` satırları geliştirme tablosuna
alınmaz ve test metriği üretilmez. Doğrulama metrikleri ile katsayılar
`reports/generated/logistic_baseline/`, eğitilmiş yerel modeller ise
`artifacts/logistic_baseline/` altında tutulur; bu iki klasör Git'e eklenmez.

Mevcut yerel veri anlık görüntüsündeki doğrulama sonuçları:

| Model | ROC-AUC | F1 | Dengeli doğruluk |
| --- | ---: | ---: | ---: |
| Eksiklik bayrağı olmadan | 0.7845 | 0.6945 | 0.7211 |
| Eksiklik bayrağıyla | 0.7850 | 0.6945 | 0.7211 |

Eksiklik bayrağının katkısı çok küçük olduğundan daha temkinli ana baseline,
bu bayrağı kullanmayan modeldir. Bu değerler doğrulama sonucudur; nihai test
performansı değildir.

### Validation karar eşiği analizi

ROC/precision-recall eğrilerini ve farklı duyarlılık hedeflerindeki karar
eşiklerini üretmek için:

```powershell
.\.venv\Scripts\python.exe scripts\analyze_validation_thresholds.py
```

Analiz yalnızca kaydedilmiş `validation` tahminlerini kabul eder; dosyada başka
bir ayrım bulunursa betik hata verir. Sabit `0.5`, en yüksek Youden J ve en az
`%80 / %90 / %95` duyarlılık çalışma noktaları
`reports/generated/logistic_baseline/threshold_analysis/` altında raporlanır.
Bir eşik validation üzerinde seçildikten sonra model ve eşik sabitlenmeli;
nihai `test` değerlendirmesi yalnızca bir kez yapılmalıdır.

Eksiklik bayrağı olmayan ana model için mevcut validation çalışma noktaları:

| Politika | Eşik | Duyarlılık | Özgüllük | Yanlış negatif | Yanlış pozitif |
| --- | ---: | ---: | ---: | ---: | ---: |
| Sabit eşik | 0.5000 | %65,12 | %79,09 | 9.272 | 5.959 |
| En az %80 duyarlılık | 0.3381 | %80,14 | %57,29 | 5.281 | 12.169 |
| En az %90 duyarlılık | 0.2838 | %90,55 | %35,45 | 2.512 | 18.394 |
| En az %95 duyarlılık | 0.2408 | %95,37 | %24,12 | 1.231 | 21.623 |

Bu tablo bir eşik önerisi değil, yanlış negatiflerle yanlış pozitifler arasındaki
ödünleşmenin validation verisindeki ölçümüdür. Kullanım senaryosunun hata
maliyetleri tanımlanmadan operasyonel eşik sabitlenmez.

### Validation alt grup hata analizi

Ana modelin sabit `0.5` eşiğindeki hatalarını varyant türü, kromozom, mutasyon
sınıfı, baz değişimi, SNV durumu ve uzunluk farkına göre incelemek için:

```powershell
.\.venv\Scripts\python.exe scripts\analyze_validation_errors.py
```

Alt grup metrikleri, en güvenli görünüp yanlış çıkan örnekler ve hata oranı
grafikleri `reports/generated/logistic_baseline/error_analysis/` altında
üretilir. Tabloda tüm gruplar korunur; 100'den az validation kaydı olan gruplar
kararsız oran uyarısıyla işaretlenir ve varsayılan grafiklere alınmaz. Betik
yalnızca validation tahmin kimliklerinin özelliklerini belleğe alır.

Mevcut validation analizindeki başlıca bulgular:

- Tek nükleotid varyantları 36.224 kayıtla en büyük gruptur; bu gruptaki
  patojenik kayıtların `%77,49`u `0.5` eşiğinde yanlış negatif olmuştur.
- Toplam 9.272 yanlış negatifin 9.152'si tek nükleotid varyantıdır.
- Transition grubundaki 6.351 patojenik kaydın tamamı yanlış negatif olmuştur.
- Deletion ve duplication gruplarında duyarlılık `%100` iken özgüllük `%0`dır;
  model bu gruplardaki bütün kayıtları patojenik sınıfa atamaktadır.
- Bu davranış, ilk baseline'ın kaba varyant kategorilerine fazla yaslandığını ve
  aynı varyant türü içindeki biyolojik ayrımı yeterince öğrenemediğini gösterir.

Bu nedenle sonraki model geliştirme adımı yalnızca genel skoru yükseltmeye
değil, özellikle SNV'lerde sınıf ayrımını güçlendirecek biyolojik özelliklere ve
özellik etkileşimlerine odaklanmalıdır.

### Gen-grubu çakışmasız veri ayrımı

Aynı genin veya ortak gen taşıyan çok-genli kayıtların farklı veri bölümlerine
düşmesini engelleyen daha katı değerlendirme ayrımını üretmek için:

```powershell
.\.venv\Scripts\python.exe scripts\build_gene_disjoint_split.py
```

Betik noktalı virgül veya virgülle ayrılmış gen sembollerini bağlı bileşenler
olarak birleştirir. Yirmi sınıf-korumalı gen katının 14'ü train, 3'ü validation,
3'ü test olarak atanır. Çıktı
`data/processed/baseline_gene_disjoint_split.tsv`, sayımlar ve bütünlük
kontrolleri `reports/generated/gene_disjoint_split/` altında tutulur. Bu ayrım,
varyant düzeyindeki ilk ayrıma göre daha zordur ve modelin eğitimde hiç görmediği
gen gruplarına genellemesini ölçer.

Mevcut veri anlık görüntüsünde 11.904 bağlı gen bileşeni oluşmuştur. Ayrım
sonucunda train bölümünde 257.035, validation bölümünde 55.068 ve kilitli test
bölümünde 55.093 varyant vardır; bölümler arasındaki gen grubu çakışması sıfırdır.

Aynı lojistik baseline'ı bu ayrım üzerinde çalıştırmak için:

```powershell
.\.venv\Scripts\python.exe scripts\train_logistic_baseline.py `
  --split data\processed\baseline_gene_disjoint_split.tsv `
  --report-dir reports\generated\logistic_baseline_gene_disjoint `
  --artifact-dir artifacts\logistic_baseline_gene_disjoint
```

Eksiklik bayrağı olmayan ana modelin iki validation protokolündeki karşılaştırması:

| Validation protokolü | ROC-AUC | F1 | Dengeli doğruluk |
| --- | ---: | ---: | ---: |
| Varyant düzeyinde ayrım | 0.7845 | 0.6945 | 0.7211 |
| Görülmemiş gen bileşenleri | 0.7630 | 0.6678 | 0.7018 |

Gene-disjoint validation'da daha düşük performans gözlenmiştir. Bu bulgu daha
zor yeni-gen genellemesiyle uyumludur; ancak iki ayrımdaki gen ve grup
dağılımları da değiştiği için düşüş tek bir nedene bağlanamaz. Bundan sonraki
model karşılaştırmalarında gene-disjoint validation daha katı ana geliştirme
ölçütü olarak kullanılmalıdır; bu aşamada test model seçimi için açılmamıştır.

### HGVS sonuç özellikleri

ClinVar `Name` alanındaki ham gen ve transkript metnini modele vermeden yalnızca
HGVS sonuç türlerini çıkarmak için:

```powershell
.\.venv\Scripts\python.exe scripts\build_hgvs_features.py
```

Çıktı `data/processed/clinvar_hgvs_features.tsv.gz` konumuna yazılır. Eklenen
özellikler HGVS dizi türü, DNA değişim türü, kodlayan bölge sınıfı, intronik
splice uzaklığı sınıfı, protein etkisi ve protein açıklaması varlığıdır. Ham
`Name`, gen sembolü ve transkript kimliği çıktıya alınmaz. Özellik dağılımları
`reports/generated/hgvs_features/` altında yalnızca gene-disjoint `train` ve
`validation` hedeflerine göre raporlanır; kilitli test hedef dağılımı rapora
alınmaz.

HGVS özellikli lojistik modeli gene-disjoint ayrım üzerinde eğitmek için:

```powershell
.\.venv\Scripts\python.exe scripts\train_logistic_baseline.py `
  --features data\processed\clinvar_hgvs_features.tsv.gz `
  --split data\processed\baseline_gene_disjoint_split.tsv `
  --feature-profile hgvs `
  --report-dir reports\generated\logistic_hgvs_gene_disjoint `
  --artifact-dir artifacts\logistic_hgvs_gene_disjoint
```

Eğitim kodunun varsayılan `baseline` profili korunur. `hgvs` profili seçilmeden
HGVS sütunları okunmaz; böylece iki özellik seti aynı model ve veri ayrımıyla
adil biçimde karşılaştırılabilir.

Gene-disjoint validation üzerindeki ana model karşılaştırması:

| Özellik seti | ROC-AUC | Average precision | F1 | Dengeli doğruluk |
| --- | ---: | ---: | ---: | ---: |
| Yapısal baseline | 0.7630 | 0.7516 | 0.6678 | 0.7018 |
| Yapısal + HGVS biyolojik sonuçları | 0.9825 | 0.9810 | 0.9196 | 0.9237 |
| Biyolojik sonuçlar + açık gösterim alanları | 0.9824 | 0.9810 | 0.9194 | 0.9235 |

Sabit `0.5` eşiğinde genel yanlış negatif sayısı 10.238'den 2.919'a, yanlış
pozitif sayısı 6.015'ten 1.231'e düşmüştür. SNV alt grubunda ROC-AUC
`0.6423 → 0.9714`, duyarlılık `%18,23 → %79,71` olmuş; yanlış negatif sayısı
10.077'den 2.500'e gerilemiştir.

Bu güçlü artış biyolojik olarak beklenen sonuç sınıflarından gelir: örneğin
stop-gained, frameshift ve splice çekirdeği etkileri patojeniklik yönünde;
synonymous ve deep-intronic sınıfları benign yönünde ağırlık almıştır. Bununla
birlikte sonuç hâlâ ClinVar validation verisine aittir. ClinVar'ın anotasyon ve
seçim yanlılıkları dış veri performansını yükseltebilir veya düşürebilir; bu
model klinik kullanıma hazır değildir. Bu validation analizi yapıldığı aşamada
gene-disjoint test henüz açılmamıştı; güncel değerlendirme durumu belgenin
başındaki "Mevcut durum" bölümündedir.

### HGVS ayrıştırıcı denetimi ve özellik ablasyonu

Ham `Name` alanındaki ayrıştırma örneklerini ve kaydedilmiş özelliklerin yeniden
üretilebilirliğini denetlemek için:

```powershell
.\.venv\Scripts\python.exe scripts\audit_hgvs_parser.py
```

Betik yalnızca gene-disjoint `train` ve `validation` bölümlerini ayrıştırır;
`test` satırları dosya hizası kontrolünden sonra atlanır. Mevcut veride 312.103
geliştirme kaydının altı HGVS alanında kaydedilmiş tabloyla uyuşmazlık yoktur.
Missense, synonymous, frameshift, stop, protein anotasyonsuz ve belirsiz
sınıflardan deterministik seçilen 60 ham örneğin basit bağımsız sözdizimi işaret
kontrolü de 60/60 tutarlıdır. Bu kontrol klinik veya transkript-duyarlı uzman
doğrulamasının yerine geçmez.

Bu denetim sırasında `Xaa` ile bilinmeyen alternatif aminoasit belirten üç
train kaydının genel üç-harfli kalıp tarafından missense sayıldığı görüldü.
Ayrıştırıcı düzeltildi ve bu kayıtlar `uncertain_or_absent` sınıfına taşındı;
özellik tabloları ve geliştirme karşılaştırmaları yeni veri hash'iyle yeniden
üretildi.

Üç özellik profilini aynı lojistik model ve aynı gene-disjoint ayrım üzerinde
karşılaştırmak için:

```powershell
.\.venv\Scripts\python.exe scripts\run_hgvs_ablation.py
```

Ablasyon, performans artışının HGVS biyolojik sonuç sınıfları eklendiğinde
oluştuğunu gösterir: yapısal profile göre ROC-AUC farkı `+0.2194`, AP farkı
`+0.2294` olmuştur. `hgvs_sequence_type` ve `has_protein_annotation` açık
gösterim alanlarının ayrıca eklenmesi ROC-AUC'yi `-0.000123`, AP'yi `-0.000063`
değiştirmiştir; pratik ek katkı görülmemiştir. Bununla birlikte biyolojik
profildeki `protein_effect=no_protein_annotation` kategorisi anotasyon yokluğunu
dolaylı taşıyabildiği için bu karşılaştırma bütün anotasyon-varlığı sinyalini
tamamen ortadan kaldıran bir deney değildir.

Tam modelin validation alt grup analizi genel skorun önemli bir zayıflığı
gizlediğini göstermiştir: 7.075 missense kaydında ROC-AUC `0.5696`, duyarlılık
`%31,49`dır. Frameshift ve stop gruplarında pozitif sınıf oranı sırasıyla
`%99,40` ve `%99,69` olduğundan çok yüksek AP/duyarlılık değerleri güçlü sınıf
dengesizliğiyle birlikte yorumlanmalıdır. Ayrıntılı, sürüm kontrollü sonuç özeti
`docs/hgvs_parser_ablation_result.md` dosyasındadır; yerel ham örnek ve tablo
çıktıları `reports/generated/hgvs_parser_audit/` ile
`reports/generated/hgvs_ablation/` altındadır.

### Aynı özelliklerle ağaç modeli karşılaştırması

Tam HGVS özellik setinde doğrusal olmayan etkileşimlerin katkısını tek bir sabit
ağaç baseline'ıyla ölçmek için:

```powershell
.\.venv\Scripts\python.exe scripts\compare_tree_model.py
```

Betik referans lojistik regresyon ile `HistGradientBoostingClassifier` modelini
aynı gene-disjoint train/validation ayrımı, aynı özellikler, `seed=2027` ve
sabit `0.5` eşikle karşılaştırır. Ağaç ayarları tek koşudan önce sabitlenmiştir;
validation üzerinde hiperparametre taraması yapılmaz ve test bölümü yüklenmez.

| Model | ROC-AUC | AP | Recall | Özgüllük | F1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Lojistik regresyon | 0.9824 | 0.9810 | 0.8902 | 0.9568 | 0.9194 |
| HistGradientBoosting | 0.9839 | 0.9821 | 0.8954 | 0.9558 | 0.9218 |

Ağaç modeli genel ROC-AUC'de `+0.00155`, AP'de `+0.00114` ve recall'da
`+0.00523` artış sağlamış; yanlış negatif sayısını 2.919'dan 2.780'e düşürürken
yanlış pozitif sayısını 1.231'den 1.259'a çıkarmıştır. Ancak missense ROC-AUC
`0.5696 → 0.5555`, missense recall `0.3149 → 0.2893` olmuştur. Synonymous
pozitifleri yine yakalanmamış; frameshift ve stop-gained/changed gruplarında
özgüllük yine sıfır kalmıştır. Bu nedenle küçük genel artış, ana modelleme
sorununun çözüldüğü anlamına gelmez. Ayrıntılı sonuçlar
`docs/tree_model_comparison_result.md` dosyasındadır.

### Hedefli missense biyolojik özellikleri

Standart missense protein değişimlerinden gen, transkript ve protein konumunu
modele taşımadan aminoasit çifti, Grantham uzaklığı, BLOSUM62 skoru, hidropati
farkı ve fizikokimyasal sınıf geçişleri üretmek için:

```powershell
.\.venv\Scripts\python.exe scripts\build_missense_features.py
```

47.708 missense kaydın tamamı standart aminoasit çifti olarak ayrıştırılmıştır.
Bu özelliklerin kontrollü lojistik regresyon ablasyonunu çalıştırmak için:

```powershell
.\.venv\Scripts\python.exe scripts\run_missense_feature_ablation.py
```

| Profil | Genel ROC-AUC | Genel F1 | Missense ROC-AUC | Missense AP | Missense recall |
| --- | ---: | ---: | ---: | ---: | ---: |
| HGVS referans | 0.9824 | 0.9194 | 0.5696 | 0.5078 | 0.3149 |
| + missense biyolojik özellikler | 0.9859 | 0.9381 | 0.7513 | 0.6878 | 0.5704 |

Kazanç missense grubunda yoğunlaşmıştır: ROC-AUC `+0.1817`, AP `+0.1800`,
recall `+0.2554` ve F1 `+0.2317`. Bununla birlikte bu hâlâ tek sabit validation
ölçümüdür; aminoasit değişim kategorileri ClinVar seçim örüntülerini de
öğrenebilir. Ana modele geçmeden önce train-only gene-disjoint çapraz doğrulama
ve daha sonra bağımsız dış değerlendirme gerekir. Ayrıntılar
`docs/missense_feature_ablation_result.md` dosyasındadır.

### Train-only gene-disjoint çapraz doğrulama

HGVS modelinin tek bir gen ayrımına özgü olup olmadığını, sabit validation ve
test bölümlerine dokunmadan ölçmek için:

```powershell
.\.venv\Scripts\python.exe scripts\cross_validate_gene_disjoint.py
```

Betik yalnızca gene-disjoint `train` bölümünde 5 katlı
`StratifiedGroupKFold` uygular. Her gen bileşeni tek bir holdout fold'unda kalır
ve her train varyantı tam bir out-of-fold tahmin alır. Fold metrikleri, birleşik
out-of-fold tahminleri, katsayı kararlılığı ve görsel özet
`reports/generated/hgvs_gene_disjoint_cv/` altında üretilir. Sabit validation
ve test bu analizde kullanılmaz.

Gerçek train koşusunda 257.035 varyant ve 8.132 gen bileşeni kullanılmıştır.
Beş fold için ROC-AUC ortalaması `0.9828 ± 0.0041` (örneklem standart sapması),
aralığı `0.9781–0.9867` olmuştur. Birleştirilmiş out-of-fold sonuçları:

| Metrik (sabit eşik = 0.5) | OOF değer |
| --- | ---: |
| ROC-AUC | 0.9828 |
| Average precision | 0.9822 |
| F1 | 0.9216 |
| Dengeli doğruluk | 0.9253 |
| Duyarlılık | 0.8978 |
| Özgüllük | 0.9529 |

Her train varyantı tam bir OOF tahmin almış, hiçbir gen bileşeni birden fazla
holdout fold'una dağılmamıştır. Bu sonuçların sabit gene-disjoint validation
ROC-AUC değeri olan `0.9824` ile yakın olması, model performansının seçilmiş tek
bir train/validation ayrımına güçlü biçimde bağlı olmadığına dair olumlu bir
kararlılık kontrolüdür. Bu analiz dış veri doğrulamasının yerine geçmez ve
bu train-only analiz çalıştırıldığında test bölümü henüz açılmamıştı. Güncel
değerlendirme durumu belgenin başındaki "Mevcut durum" bölümündedir.

### Missense özelliği için train-only OOF kararlılık kontrolü

Missense biyolojik özelliklerinin validation kazanımının tek bir ayrımın tesadüfü
olup olmadığını sınamak için, iki profil aynı train satırlarında ve aynı 5-fold
`StratifiedGroupKFold` atamalarıyla tekrar karşılaştırıldı:

```powershell
.\.venv\Scripts\python.exe scripts\run_missense_gene_disjoint_cv.py
```

Sabit validation ve kilitli test bu koşuda yüklenmez. Birleştirilmiş OOF sonuçları
aynı yönü doğrulamıştır:

| Grup | HGVS ROC-AUC | + missense ROC-AUC | HGVS AP | + missense AP | HGVS recall | + missense recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Tüm train OOF | 0.9829 | 0.9862 | 0.9823 | 0.9857 | 0.8959 | 0.9262 |
| Missense train OOF | 0.5503 | 0.7358 | 0.4939 | 0.6811 | 0.3200 | 0.5691 |

Missense alt grubundaki farklar ROC-AUC için `+0.1855`, AP için `+0.1871` ve
recall için `+0.2491` olmuştur. Aynı satırların fold ataması iki profil arasında
VariationID düzeyinde doğrulanmıştır. Bu, gözlemin yalnızca tek sabit validation
ayrımına bağlı olmadığını destekleyen bir iç geliştirme kanıtıdır; bağımsız dış
veri ya da klinik geçerlilik kanıtı değildir. Ayrıntılar
`docs/missense_gene_disjoint_cv_result.md` dosyasındadır.

### Dış doğrulama adayı modeli dondurma

Profil seçimi tamamlandıktan sonra aday modeli yalnızca `train + validation`
üzerinde yeniden eğitip, dış kohort etiketleri görülmeden dondurmak için:

```powershell
.\.venv\Scripts\python.exe scripts\freeze_external_candidate.py
```

Bu işlem `hgvs_missense` profilini, sabit `0.5` karar eşiğini ve artifact
SHA-256 değerini yerel model kartına yazar. Kilitli iç test kullanılmaz. Dış
kohortun kabul kriterleri ve değerlendirme kuralı
`docs/external_validation_protocol.md` dosyasında tanımlıdır.

Dış ham kohort, önce geliştirme verisinden çakışmasız hâle getirilmelidir:

```powershell
.\.venv\Scripts\python.exe scripts\admit_external_cohort.py `
  --external data\external\candidate_binary.tsv.gz `
  --cohort-kind independent `
  --source-version "sağlayıcının-sürümü" `
  --source-url "https://kaynak.example.org" `
  --reference-config data\reference\grch38.json
```

Bu adım `Assembly=GRCh38`, referans alleli, `VariationID`, `AlleleID` ve referans
üzerinde sola hizalanmış koordinat-allel çakışmalarını denetler. Desteklenmeyen
veya eksik alleller dışlanır; dış kohortun kendi içindeki eşdeğer tekrarlar
işlemi durdurur. Referans FASTA ve FAI hash'leri `grch38.json` içinde sabitlenir;
örnek config şeması `docs/external_validation_protocol.md` içindedir.

Kabul edilmiş dış ikili kohorttan aynı özellik şemasını üretmek için:

```powershell
.\.venv\Scripts\python.exe scripts\build_external_features.py `
  --input data\external\admitted_binary.tsv.gz `
  --output data\external\admitted_missense_features.tsv.gz
```

Ardından dondurulmuş artifact ile tek-seferlik ölçüm için:

```powershell
.\.venv\Scripts\python.exe scripts\evaluate_external_candidate.py
```

Değerlendirme etiketleri okumadan önce kalıcı `EVALUATION_STARTED.json` ve aday
model klasöründe kohort hash'ine bağlı bir registry kaydı oluşturur. Yarım kalan
koşu farklı rapor klasörüyle de yeniden başlatılmaz. Rapor tüm kohort, missense
ve görülmemiş genler için metrikleri, karar hata sayılarını ve en az 1.000 tekrarlı
bootstrap güven aralıklarını içerir. İki sınıfı bulunmayan alt gruplar açıkça
işaretlenir. Yapay verili uçtan uca denetim:

```powershell
.\.venv\Scripts\python.exe -m unittest tests.test_external_pipeline
```

Açık kaynak adaylarının ön incelemesi: `docs/external_data_source_review.md`.

### ClinGen eRepo kimlik ve uygunluk taraması

İlk açık kaynak adayı için resmî arayüzün kullandığı API'den yerel anlık görüntü
almak ve eski 367.196 varyantla kimlik örtüşmesini ölçmek için:

```powershell
.\.venv\Scripts\python.exe scripts\download_clingen_erepo.py `
  --output-dir data\external\clingen_erepo_2026-10-26
.\.venv\Scripts\python.exe scripts\screen_clingen_erepo.py `
  --snapshot-dir data\external\clingen_erepo_2026-10-26
```

İndirme klasörü doluysa üzerine yazılmaz. Sayfaların ve kullanım şartlarının
hash'leri saklanır; UUID tekilliği ve indirme öncesi/sonrası toplam sayılar
karşılaştırılır. Bu, canlı API'nin indirildiği anda sabitlenmiş yerel kopyasıdır;
sağlayıcının resmî sürümlü veri arşivi değildir.

Tarama `reports/generated/clingen_erepo_screen/` altında tüm etiket sayılarını,
çelişen yorumları ve eski VariationID ile eşleşmeyen adayları yazar. Hastalık,
uzman paneli ve yorum sürümleri korunur. `Likely` etiketleri kesin ikili hedefe
katılmaz. Eşleşmeyen kimlikler, allel/koordinat ve etiket-kökeni denetimleri
tamamlanana kadar yalnızca aday üst sınırıdır. Bu aşamada model yüklenmez.

26 Ekim 2026 taramasında 13.564 yorum kaydından, tek ClinVar kimliği bulunan
4.610 kesin ikili varyant belirlendi. Eski veriyle 4.321 kimlik çakışması ve 7
uyuşmayan sınıflandırma çıkarılınca 282 aday (216 Pathogenic, 66 Benign) kaldı.
Bunların 180'i missense olarak ayrışır. Henüz allel/koordinat kabulü yapılmadı.
Deneyin sabitlenmiş özeti: `docs/clingen_erepo_screen_result.md`.

### OOF olasılık kalibrasyonu

Gene-disjoint çapraz doğrulamadaki olasılıkların güvenilirliğini incelemek için:

```powershell
.\.venv\Scripts\python.exe scripts\analyze_oof_calibration.py
```

Betik yalnızca train out-of-fold tahminlerini kabul eder. Eşit-frekans ve sabit
genişlikli olasılık dilimlerini, fold bazında kalibrasyon ölçülerini ve
güvenilirlik grafiğini
`reports/generated/hgvs_gene_disjoint_cv/calibration/` altında üretir. Bu aşama
bir tanısal analizdir; modele yeniden kalibrasyon uygulamaz ve sabit validation
veya test verisini kullanmaz.

257.035 OOF tahmini üzerindeki sonuçlar:

| Kalibrasyon ölçüsü | Değer |
| --- | ---: |
| Brier skoru | 0.0466 |
| Log loss | 0.1525 |
| ECE (10 eşit-frekans dilimi) | 0.0050 |
| MCE (10 eşit-frekans dilimi) | 0.0220 |
| Kalibrasyon intercept'i (ideal 0) | 0.0213 |
| Kalibrasyon slope'u (ideal 1) | 0.9662 |

Genel güvenilirlik eğrisi ideale yakındır. Bununla birlikte sabit genişlikli
orta olasılık dilimlerinde yerel sapmalar vardır: `0.50–0.60` bandında tahmin
gözlenenden yaklaşık 5,4 yüzde puanı yüksek, `0.70–0.80` bandında yaklaşık 6,0
yüzde puanı düşüktür. Bu nedenle genel bir yeniden kalibrasyon şimdilik
uygulanmamış; model olasılıklarının bireysel klinik risk olarak değil, öncelikle
sıralama ve araştırma amaçlı skor olarak ele alınması korunmuştur.

### OOF belirsiz bölge ve hata analizi

Sabit `0.5` karar eşiğinin çevresindeki `0.20–0.80` olasılık bölgesini HGVS ve
varyant alt gruplarında incelemek için:

```powershell
.\.venv\Scripts\python.exe scripts\analyze_oof_uncertainty.py
```

Betik train OOF tahminlerini HGVS özellikleriyle VariationID üzerinden
eşleştirir. Düşük, belirsiz ve yüksek olasılık bölgelerinin hata yükünü;
belirsiz bölgedeki alt grup hata oranlarını ve belirsiz bölge dışında kalan
yüksek güvenli hataları raporlar. Validation ve test verisi kullanılmaz.

Gerçek OOF analizinde `0.20–0.80` bandı 39.924 varyantı (tüm train OOF
verisinin `%15,53`'ü) kapsamış ve 18.943 hatanın 16.223'ünü (`%85,64`)
yakalamıştır. Bölgelere göre hata yükü:

| Olasılık bölgesi | Varyant | Hata | Hata oranı |
| --- | ---: | ---: | ---: |
| `<0.20` | 112.155 | 1.946 yanlış negatif | %1,74 |
| `0.20–0.80` | 39.924 | 16.223 | %40,63 |
| `>0.80` | 104.956 | 774 yanlış pozitif | %0,74 |

Belirsiz bölge ağırlıklı olarak SNV/substitution ve missense varyantlarından
oluşur; missense grubunda hata oranı `%44,63`'tür. Belirsiz bölge dışındaki
1.946 yanlış negatifin 1.725'inde protein anotasyonu yoktur ve 1.357'si
introniktir. 774 yüksek olasılıklı yanlış pozitifte frameshift (216) ve
stop-gained/changed (143) sonuçları öne çıkar. Bu durum modelin genel sonuç
sınıflarını güçlü öğrendiğini, fakat gen/transkript bağlamına özgü istisnaları
tek başına çözemediğini gösterir.

`0.20–0.80` bandı manuel inceleme kuyruğunu önceliklendirmek için yararlıdır;
ancak düşük olasılık bölgesinde 1.946 patojenik varyant kaldığı için bandın
dışındaki tahminler otomatik klinik karar veya “güvenli varyant” etiketi olarak
kullanılamaz.

### Test öncesi protokolü dondurma

Ana model, validation çalışma eşiği, veri/artifact hash'leri ve test başarı
kapılarını test etiketlerini okumadan kilitlemek için:

```powershell
.\.venv\Scripts\python.exe scripts\freeze_evaluation_protocol.py
```

İnsan tarafından belirlenen kararlar
`config/final_evaluation_protocol.json` içinde tutulur. Betik tüm geliştirme
manifestlerinin test kullanmadığını, veri hash'lerinin birbiriyle tutarlı
olduğunu ve model artifact'inin değiştirilmediğini doğrular. Ardından
`reports/generated/final_protocol/` altında makine-okunur kilit ve protokol
kartı oluşturur. Bu komut testi çalıştırmaz veya test etiketlerini okumaz.

Dondurulan protokol sürümü `1.0` için ana model
`without_missingness_flag`, HGVS özellik profili ve `0.3623531848` karar
eşiğidir. Bu eşik validation üzerinde `%95,02` duyarlılık ve `%88,82`
özgüllük sağlamıştır. Kilitli testin başarılı sayılması için aşağıdaki
koşulların tümü birlikte aranacaktır:

| Başarı kapısı | Sınır |
| --- | ---: |
| ROC-AUC | ≥ 0.97 |
| Average precision | ≥ 0.97 |
| Donmuş eşikte duyarlılık | ≥ 0.93 |
| Donmuş eşikte özgüllük | ≥ 0.85 |
| Validation'a göre ROC-AUC düşüşü | ≤ 0.02 |

Test yalnızca bir kez açılır. Sonuç görüldükten sonra model veya eşik
değiştirilmeden birincil sonuç raporlanır; başarısızlık durumunda gelecekteki
iddialar için yeni bağımsız holdout gerekir.

### İç gene-disjoint değerlendirme sonucu

Gene-disjoint test 55.093 varyant ve son modelin eğitiminde görülmemiş 1.798 gen
grubunda tek sefer değerlendirilmiş, önceden belirlenen beş proje-içi başarı
kapısının tamamı geçilmiştir:

| Test metriği | Değer |
| --- | ---: |
| ROC-AUC | 0.9872 |
| Average precision | 0.9870 |
| Duyarlılık | 0.9661 |
| Özgüllük | 0.8762 |
| F1 | 0.9206 |
| Dengeli doğruluk | 0.9212 |

Donmuş `0.3623531848` eşiğinde 902 yanlış negatif ve 3.527 yanlış pozitif
gözlenmiştir. Test ROC-AUC değeri validation değerinden yaklaşık 0.0048 daha
yüksektir.

İlk koşudaki `chromosome` karışık veri tipi uyarısı şeffaf biçimde ayrı bir
teknik denetimle incelenmiştir. Eğitimdeki string kategorik veri sözleşmesi
uygulandığında 208 olasılık değişmiş, hiçbir sınıf kararı değişmemiş ve bütün
başarı kapıları yine geçilmiştir. Orijinal koşu korunmuştur. Ayrıntılı sonuç ve
bütünlük hash'leri `docs/final_test_result.md` dosyasındadır. İlk ve düzeltilmiş
ROC-AUC/AP değerleri ayrı raporlanır; sınıf kararlarının aynı kalması sıralama
ve kalibrasyon ölçülerinin tam olarak aynı kaldığı anlamına gelmez.

Bu başarı kapıları proje başlamadan tanımlanmış klinik veya TEKNOFEST standartları
değil, bu çalışma içinde test açılmadan önce belirlenen iç kabul ölçütleridir.
Yarışma başarısı veya klinik geçerlilik standardı olarak sunulamazlar.

### Eski/yeni ayrım geçiş denetimi

Varyant-düzeyi ilk ayrımla gene-disjoint ayrımın rollerini karşılaştırmak için:

```powershell
.\.venv\Scripts\python.exe scripts\audit_split_transition.py
```

Son gene-disjoint test varyantlarının önceki rolleri:

| Önceki varyant-düzeyi rol | Varyant | Son test içindeki oran |
| --- | ---: | ---: |
| Train | 38.630 | %70,12 |
| Validation | 8.200 | %14,88 |
| Test | 8.263 | %15,00 |

Son model doğrudan yalnızca yeni gene-disjoint train bölümünde eğitildi ve gen
grubu çakışması yoktu. Buna karşın son testin 46.830 varyantı (`%85,00`) daha
önceki deneylerde train veya validation rolünde bulundu. Önceki validation
sonuçları hata incelemesini ve özellik geliştirmeyi etkilediği için mevcut test
geliştirme süreci açısından tamamen görülmemiş sayılamaz. Sonuç bu nedenle
**iç gene-disjoint değerlendirme** olarak sunulur; bağımsız nihai iddia için
harici değerlendirme gerekir.

