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

## Veri

Ham veri büyük olduğu ve harici bir kaynaktan yeniden indirilebildiği için Git
repository'sine eklenmez. Beklenen yerel konum:

```text
data/raw/variant_summary.txt.gz
```

Kaynak: [NCBI ClinVar indirme sayfası](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/)

Doğrudan dosya: [variant_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)

## Sıradaki adım

Ham dosyayı parça parça okuyarak genom referansı, klinik etiket, köken,
değerlendirme durumu, eksik değer ve tekrar dağılımlarını raporlayan yeniden
çalıştırılabilir bir veri denetim betiği hazırlanacak.

