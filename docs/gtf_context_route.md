# Mevcut diskte VEP: GTF tabanlı çevrimdışı yol

Tam insan VEP cache kurulumu için ayrılan yaklaşık 64 GB bütçe, o kurulum
yönteminin planıydı. Ek disk veya sunucu bulunmadığı için daha küçük bir kaynak
yolu hazırlandı: **Ensembl 116 GTF + mevcut GRCh38 FASTA + sabit VEP 116.1 imajı**.

VEP'in resmî [GFF/GTF kaynak rehberi](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff)
bu dosyalarla yerel transkript modeli kurulmasını destekler. Aynı rehber,
GTF-only kullanımda protein domainleri gibi yardımcı bilgilerin eksik
kalabileceğini açıklar. Bu nedenle GTF yolu, önceki cache/domain profilinden
ayrı tutulur; eksik Pfam/InterPro değerleri sıfırla doldurulmaz.

## Kaynaklar ve disk kullanımı

GTF kaynağı:

```text
https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
```

Sıkıştırılmış kaynak **141.121.632 bayt** (141,12 MB). Başlığı GRCh38.p14
anotasyonunu belirtir. İndirme aralıkları HTTP 206, Content-Range, toplam boyut ve
ETag eşleşmesiyle kabul edilir; her parça SHA256 makbuzuyla korunur. Eksik veya
sürümü karışmış aralıklar birleştirilmez. Kaynak dosya, sıralı bgzip GTF, tabix
indeksi ve transkript SQLite verisi yerelde hash ile kilitlenir.

İlk sıralama sırasında geçici dosyalar gerektiğinden **8 GiB hazırlık bütçesi**
ayrıldı. Bu, kalıcı kaynak boyutu değildir. Kurulumdan sonraki çalışma payı
2 GiB'dir; her yeni VCF parçasından önce boş alan kontrol edilir. Tam cache
arşivi indirilmez; mevcut FASTA'nın ikinci kopyası oluşturulmaz.

Varyant kapsamı önceki hazırlıkla aynıdır: yalnızca train bölümündeki 169.265
birincil kromozom SNV'si, 170 parça. Aynı VCF dosyaları yeniden kullanılabilir.
Bu yol indelleri, yapısal varyantları veya alternatif contigleri kapsamaz.

## Komutlar

İlk kaynak indirme ve kurulum:

```powershell
.\.venv\Scripts\python.exe scripts\prepare_gtf_resources.py --fetch
```

Sonraki çağrıda `--fetch` kaldırıldığında kaynaklar ağ kullanmadan doğrulanır.
GTF sıralama, bgzip/tabix ve tüm VEP anotasyonları ağ erişimi kapalı konteynerde
yürür. Python indirmesi dışında dış servis sorgusu gerekmez.

VEP 116'nın `--offline` kısayolu cache kullanımını zorla açtığından GTF-only
komutunda bu kısayol kullanılmaz: `--cache 0 --gtf ...` ile yalnızca dosya
kaynağı seçilir, `--database` açılmaz ve ağ sınırı Docker `--network=none` ile
uygulanır. Bu davranış gerçek konteynerde sınandı; GTF yolu tam cache aramaz.

Kaynak kilidi ve ortam kontrolü:

```powershell
.\.venv\Scripts\python.exe scripts\local_vep_pipeline.py lock-resources --config config\local_vep_gtf.json
.\.venv\Scripts\python.exe scripts\local_vep_pipeline.py check --config config\local_vep_gtf.json
```

Önceki REST pilotunun 15 kaydını aynı koordinatlarla yerelde anotasyonlamak ve
özellik eşleşmesini incelemek:

```powershell
.\.venv\Scripts\python.exe scripts\compare_gtf_context_pilot.py
```

Bu, teknik karşılaştırmadır; uzman etiketlerini veya model skorlarını okumaz.
Önceden görülmüş üç tanısal örnek yeni bağımsız doğrulama olarak sunulmaz.

Train parçalarını çalıştırmak ve tamamlanan parçalardan özellik çıkarmak:

```powershell
.\.venv\Scripts\python.exe scripts\local_vep_pipeline.py run --config config\local_vep_gtf.json --max-chunks 1
.\.venv\Scripts\python.exe scripts\compare_gtf_context_pilot.py --completed-batches
```

`--max-chunks 170` kalan parçaları işler; doğrulanmış tamamlanma makbuzu olanlar
yeniden hesaplanmaz. Girdi, kaynak, kod/ayar veya çıktı değişirse sessizce
karıştırmak yerine işlem durur. Başarısız denemeler ayrı dizinde korunur.

## Özelliklerin kaynağı ve sınırı

VEP GTF'den varyantın transkript sonuçlarını üretir. SQLite tablosu aynı GTF'den
gen sembolü, transkript sürümü, MANE_Select/Ensembl_canonical işaretleri, strand,
ekzonlar ve CDS uzunluğunu taşır. Transkript seçim sırası önceki pilotla aynıdır:
gen sembolü eşleşmesi, MANE Select, canonical, alfabetik transkript kimliği.
Önceki pilotun seçimini sonuçlara bakarak kopyalayan bir kural kullanılmaz.

Aday biyolojik sütunlar:

- `biotype`, `consequence_terms`;
- `exon_rank`, `exon_total`, `location`;
- `donor_distance_bp`, `acceptor_distance_bp`, `nearest_splice_boundary_bp`;
- `protein_relative_position`.

Protein göreli konumu, VEP protein başlangıcı / GTF CDS uzunluğundan hesaplanan
protein uzunluğudur. Tam olmayan CDS işareti veya uygun uzunluk bulunmadığında
boş kalır. Splice uzaklıkları önceki strand-duyarlı hesapla üretilir; transkript
uçları splice birleşimi sayılmaz. RNA deneyi veya SpliceAI skoru üretilmez.

Pfam/InterPro sütunları denetim tablosunda boş bırakılır ve aday özellik
tablosundan çıkarılır. Gen/transkript kimlikleri denetim alanlarıdır; VariationID
yalnızca birleştirme anahtarıdır. Etiket ve klinik kanıt metinleri okunmaz.

## Yerel çıktılar

- `data/reference/ensembl_gtf_116/`: sıkıştırılmış GTF, sıralı GTF, indeks, SQLite ve kaynak makbuzları.
- `reports/generated/gtf_context_pilot/`: yerel pilot, REST ile sütun karşılaştırması, farklılıklar ve manifest.
- `reports/generated/local_vep_gtf_batches/`: ham parça anotasyonları ve tamamlanma makbuzları.
- `reports/generated/local_vep_gtf_features/`: yalnızca tamamlanan train parçalarının denetim ve özellik tabloları; manifest tüm 169.265 kaydın tamamlanıp tamamlanmadığını ayrıca belirtir.

Model ve ClinGen birincil sonuç hash'leri korunur. Bu kaynak değişikliğinin
sınıflandırma başarısına etkisi henüz ölçülmez; sonraki model deneyi train-only
gene-disjoint geliştirme düzeni ve yeni dış değerlendirme verisi gerektirir.

## 27 Ekim 2026 sonucu (Türkiye saati)

**Ek disk/sunucu kullanmadan GTF yolu çalıştırıldı.** Kaynaklar kuruldu,
kilitlendi ve ortam kontrolü `ready: true` verdi.

| Kaynak / çıktı | Ölçülen değer |
| --- | ---: |
| GTF'den indekslenen transkript | 644.255 |
| Kilitlenen beş kaynak dosyası | 733.226.628 bayt, yaklaşık 733 MB |
| İndirme parçaları ve tanılama dosyalarıyla tüm yeni kaynak dizini | 878.618.742 bayt, yaklaşık 879 MB |
| REST ile karşılaştırılan önceden görülmüş örnek | 15 |
| Yerelde tamamlanan train parçası | 1 / 170 |
| Gerçek train anotasyonu | 1.000 / 169.265 |
| İlk parçada splice geometrisi bulunan | 972 / 1.000 |
| İlk parçada protein göreli konumu bulunan | 496 / 1.000 |

733/879 MB ölçümleri mevcut FASTA ve daha önce kurulan Docker imajını içermez.
Hazırlık sırasındaki geçici sıralama dosyaları kalıcı kaynak boyutu değildir.

15 örneğin tamamında seçilen transkript ve dokuz biyolojik sütun, kaydedilmiş
REST pilotuyla eşleşti: biotype, sonuç terimleri, ekzon sıra/toplamı, konum
sınıfı, üç splice uzaklığı ve protein göreli konumu. Eksik değerler de karşılaştırmaya
katıldı; örneğin intronik protein konumlarının iki tarafta da boş olması eşleşmedir.
Pfam/InterPro bu karşılaştırmanın dışında tutuldu. Yeni VEP skorları veya model
tahminleri üretilmedi.

Bu eşleşme için etiket veya REST değerlerini kopyalayan düzeltme yapılmadı.
Biotype ve seçim işaretleri aynı GTF metadatasından, varyant sonuçları yerel VEP'ten,
ekzon geometrisi/CDS uzunluğu yerel transkript tablosundan alındı. VEP'in GTF
yolunda ayrıca bulunan UTR/start_codon satırlarını desteklemediğine dair
uyarıları korundu; modeller exon/CDS/stop_codon kayıtlarıyla kuruldu. 15 örnekteki
uyum tüm gen/transkriptlerde birebir eşitlik iddiası değildir.

İlk train parçası genom sıralı olduğundan yalnızca kromozom 1'in başlangıcını
kapsar. Exonic 530, intronic 442, seçilen transkript dışında 28 kayıt vardır.
Splice uzaklığı boş kalan 28 kayıt son gruptadır. Bunlar atlanmadı; boş özellik
olarak korundu. Kapsam oranları tüm 169.265 varyanta genellenemez.

Bir parçanın VEP çalışma süresi yaklaşık 29 saniye, ham JSONL boyutu 11,3 MB
oldu. Bu süre kaynak hash doğrulamasını kapsamaz ve sonraki parçalar için
garanti değildir. Henüz **169 parça / 168.265 kayıt** işlenmedi; özellik manifesti
`complete: false` ile bu durumu açıkça belirtir.

Ağ çağrısı ve Docker anotasyonunu yeniden çalıştırma girişimleri hata verecek
şekilde engellenerek pilot/özellik çıkarımı tekrarlandı. REST karşılaştırması
yine farklılık göstermedi; train denetim ve özellik TSV'lerinin hash'leri birebir
aynı kaldı. Dondurulmuş model ve ClinGen sonuç hash'leri korundu.

**17 ilgili test geçti:** GTF metadatası, eksik CDS, kimlik/sürüm/strand,
biotype kaynağı, domain eksikliği, HTTP aralık bütünlüğü ve yeniden kullanım,
kısmi özellik ihracı ile önceki parçalı hattın 10 kontrolü birlikte doğrulandı.

### Bu koşunun bütünlük kayıtları

| Dosya | SHA256 |
| --- | --- |
| GTF `resources.lock.json` | `e1e29d4dc85a21326b71133753093fb860163cba387fc09160cb6439227bd1ed` |
| Pilot `manifest.json` | `8abd67c771a44147f8577be9d6d135d908173e9a9ba728252636fb0247173494` |
| Train özellik `manifest.json` | `3aea3f8a3cb810019b9b0b584157cf04b94516c513ae583dd6b260a86d90e1a6` |
| `train_context_features.tsv` | `09b78dd0a24baa60feb9b3a0054bb05a864a60b97f0f49351b0b096d1ef8b703` |

Sonraki adım kalan parçaları aynı kaynak kilidiyle işlemek, kapsam/eksikliği
tüm train SNV'leri için raporlamak ve bundan sonra aynı kayıt/fold düzeninde
özellik ekleme deneyine geçmektir. İlk 1.000 kayıt üzerinden model eğitimi veya
genel başarı iddiası yapılmadı.
