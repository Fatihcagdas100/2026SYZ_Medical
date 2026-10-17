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
- Henüz veri temizleme veya model eğitimi yapılmadı.

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

