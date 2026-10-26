# Dış doğrulama hattı teknik denetimi

Tarih: 26 Ekim 2026

## Tamamlanan düzeltmeler

- Kabul için `Assembly=GRCh38` ve SHA-256 ile sabitlenmiş FASTA/FAI zorunlu.
- Basit iki allelli varyantlar referans kontrolünden, sola hizalama ve ortak
  ekleri kaldırma işleminden geçer. Farklı yazılmış eşdeğer indeller aynı
  çakışma anahtarına iner; yanlış REF ve yanlış assembly dışlanır.
- Dış kohort içindeki eşdeğer tekrarlar işlemi durdurur.
- Değerlendirme etiketleri okumadan önce hem çıktı klasörüne başlangıç kaydı
  hem de aday klasörüne ham-kohort hash'ine bağlı tekil registry kaydı yazar.
  Yarım kalan koşular farklı rapor klasörüyle de yeniden başlatılamaz.
- Hash zinciri eksik veya değiştirilmiş dosyaları reddeder; kabul edilen ID ve
  hedeflerin özellik tablosuna aynen geçtiği ayrıca doğrulanır.
- Missense ve görülmemiş genler için metrikler, karar hata sayıları ve en az
  1.000 tekrarlı güven aralıkları eklenmiştir. İki sınıfı olmayan alt gruplar
  açık bir 'hesaplanamadı' durumu ve sayımlarıyla raporlanır.
- Dondurulmuş aday dosyalarının üzerine yazma engellenmiştir.

## Denetim kanıtı

`tests/test_external_pipeline.py`, küçük yapay FASTA, geliştirme ve dış kohortla
kabul → özellik üretimi → yapay model dondurma → değerlendirme hattını gerçek
fonksiyonları çağırarak sınar. Başarılı rapor, değiştirilmiş model/özellik dosyası,
eksik hash, kesilen koşu, ikinci koşu ve tek sınıflı alt grup senaryoları kapsanır.

`tests/test_variant_normalization.py`, tekrarlı bölgedeki eşdeğer indelleri,
FAI ile satır sınırı okumasını, referans hash bozulmasını, normalizasyonun
haplotipi korumasını ve ikinci uygulamada değişmemesini sınar.

Yapılan ölçümler yalnızca yazılım denetimidir; yapay veri skorları gerçek dış
performans olarak kullanılmaz. Gerçek aday model yeniden eğitilmedi; iç testte
yeni performans ölçümü yapılmadı.

## Gerçek değerlendirme öncesinde gereken girdiler

GRCh38 FASTA/FAI ve accession denetimli hash config'i 26 Ekim 2026'da hazırlandı.
ClinGen'den 264 varyant uzman-kürasyonlu ek değerlendirme kapsamında teknik
kabulden geçti; sonuç `docs/clingen_grch38_admission_result.md` içindedir.
Kaynak-bağımsız klinik kohort henüz doğrulanmadı.
Normalizasyon küçük iki allelli A/C/G/T varyantları kapsar;
sembolik varyantlar, liftover ve karmaşık haplotip çözümü kapsam dışıdır.
