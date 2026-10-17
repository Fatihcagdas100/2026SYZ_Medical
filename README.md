# 2026SYZ

TEKNOFEST Sağlıkta Yapay Zeka 2027 hazırlığı kapsamında genetik varyantların
patojenik/benign sınıflandırılması üzerine yürütülen çalışma.

## Mevcut durum

- Projeye özel Python 3.13 sanal ortamı `.venv` altında oluşturuldu.
- ClinVar `variant_summary.txt.gz` dosyasının başlığı ve üç örnek kaydı,
  dosyanın tamamı belleğe alınmadan incelendi.
- Dosyada 43 sütun bulundu.
- Aynı varyantın GRCh37 ve GRCh38 için ayrı satırlarda bulunabildiği görüldü.
- Birleşik klinik etiketlerin (ör. `Pathogenic/Likely pathogenic`) ayrıca ele
  alınması gerektiği belirlendi.
- Katı `GRCh38 + germline + tam Pathogenic/Benign` politikasıyla 367.196
  benzersiz varyanttan oluşan ilk veri seti hazırlandı.
- Sızıntı riski yüksek kimlik, gen ve konum alanlarını kullanmayan yapısal
  özellikler üretildi ve sabit `%70 / %15 / %15` veri ayrımı oluşturuldu.
- İlk lojistik regresyon baseline'ı yalnızca `train` üzerinde eğitildi ve
  `validation` üzerinde değerlendirildi; `test` bölümü kilitli tutuluyor.

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

Kimlik, gen ve genom konumu gibi ezberleme riski yüksek alanları dışarıda
tutarak ilk yapısal özellik tablosunu üretmek için:

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

Görülmemiş genlerdeki düşüş, yalnızca varyant-düzeyi ayrımın yeni genlere
genelleme performansını iyimser gösterebildiğini belirtir. Bundan sonraki model
karşılaştırmalarında gene-disjoint validation daha katı ana geliştirme ölçütü
olarak kullanılmalıdır; test bölümü model seçimi için açılmaz.

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
| Yapısal + HGVS sonuçları | 0.9824 | 0.9810 | 0.9195 | 0.9236 |

Sabit `0.5` eşiğinde genel yanlış negatif sayısı 10.238'den 2.902'ye, yanlış
pozitif sayısı 6.015'ten 1.241'e düşmüştür. SNV alt grubunda ROC-AUC
`0.6423 → 0.9714`, duyarlılık `%18,23 → %79,71` olmuş; yanlış negatif sayısı
10.077'den 2.500'e gerilemiştir.

Bu güçlü artış biyolojik olarak beklenen sonuç sınıflarından gelir: örneğin
stop-gained, frameshift ve splice çekirdeği etkileri patojeniklik yönünde;
synonymous ve deep-intronic sınıfları benign yönünde ağırlık almıştır. Bununla
birlikte sonuç hâlâ ClinVar validation verisine aittir. ClinVar'ın anotasyon ve
seçim yanlılıkları dış veri performansını yükseltebilir veya düşürebilir; bu
model klinik kullanıma hazır değildir ve kilitli test henüz açılmamıştır.

