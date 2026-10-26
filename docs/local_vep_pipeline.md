# Yerel VEP: parçalı eğitim verisi ve çalışma hattı

26 Ekim 2026: **169.265 train SNV'si 170 VCF parçasına hazırlandı.** VEP 116.1
konteyneri indirildi ve ağ kapalıyken sürümü/parametreleri kontrol edildi.
İnsan GRCh38 önbelleği henüz kurulu olmadığı için gerçek veride yerel toplu
anotasyon bu cache yöntemiyle çalıştırılmadı. Ek disk/sunucu bulunmadığı için
sonraki adımda [GTF tabanlı yol](gtf_context_route.md) kuruldu ve ilk 1.000 train
kaydı yerelde anotasyonlandı. Aşağıdaki disk engeli tam cache yöntemine aittir.

## Hazırlanan verinin kapsamı

| Kalem | Kayıt |
| --- | ---: |
| Gene-disjoint train | 257.035 |
| Referansla doğrulanmış, birincil kromozom, iki allelli SNV | 169.265 |
| Tek baz A/C/G/T değişimi kapsamına girmeyen | 87.098 |
| Kromozomu 1–22, X, Y dışında olan | 671 |
| Birincil kromozom referans erişimiyle eşleştirilemeyen | 1 |
| Toplam kapsam dışı | 87.770 |

169 parça 1.000'er varyant, son parça 265 varyant içerir. Sıra kromozom,
konum ve VariationID ile belirlenir. Kapsam dışı kayıtlar veri setinden silinmez;
bu çalışma hattının ayrı raporunda tutulur. SNV altkümesi için elde edilecek
sonuçlar tüm varyant türlerine genellenemez.

Hazırlık kaynak/split dosyalarının hash'lerini ve birebir ID eşleşmesini
doğrular. Yalnızca `train` seçilir; `target` sütunları okunmaz. Validation, test
ve önceki üç ClinGen tanısal örneği bu pakete katılmaz. `REF` sabit GRCh38
FASTA'dan denetlenir. Referans, assembly veya accession/kromozom uyuşmazlığı
işlemi durdurur. Eşlenemeyen alternatif contig erişimi açık kapsam dışı nedenidir.

Yerel hazırlık dizini `data/processed/local_vep_train_snv/`:

- `chunk_000001.vcf` … `chunk_000170.vcf`: etiketsiz VCF parçaları.
- `variant_map.tsv`: VariationID, gen sembolü ve doğrulanan allel/koordinatlar; denetim/birleştirme tablosu.
- `excluded_train_rows.tsv`: kapsam dışı VariationID ve nedeni.
- `chromosome_synonyms.tsv`: VCF/cache kromozom adları ile UCSC FASTA `chr` adları arasındaki eşleme.
- `manifest.json`: kaynak, referans, hazırlık kuralı, parça ve tablo hash'leri; tüm paydalar.

Gen sembolü ve tam koordinatlar bu aşamada model özelliği değildir. Parça
VCF'lerine etiket, hastalık sınıfı veya uzman kanıt metni yazılmaz.

Bu hazırlık koşusunun bütünlük kayıtları:

| Dosya | SHA256 |
| --- | --- |
| Hazırlık `manifest.json` | `5dd74eb4aaf10c3108c09bdec6b705a98cd33c11eeff2605c96ae196916fbcc7` |
| `chunk_000001.vcf` | `35480f16d54236efd945d5b672e4107e05a3cfb401a7a8f75fcda92efd68c8cd` |
| `chunk_000170.vcf` | `5251d86feeefa7b62245c784be3e00fa150704651c0b769fe5dd44de47613e5f` |
| `variant_map.tsv` | `4b8e7861f18cf76fd81655e5e624b0edfbe3ff7cfbaae668bdfca3ec512e2bb0` |

## Sürüm ve kaynak sınırı

İmaj etiketi bilgi amacıyla `release_116.1`; çalışma kimliği değişmez digest'tir:

```text
ensemblorg/ensembl-vep@sha256:d4f7ab1dad8c350b3ffbda026c1a3086310a6311b53fb90842eaf71c23801def
```

Platform `linux/amd64`; konteynerde gerçek `vep --help` çıktısı VEP `116.1` ve
Ensembl API bileşenleri `116` olarak kaydedildi. Resmî [VEP kurulum rehberi](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html)
konteyner kullanımını anlatır; [cache rehberi](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
araç ve cache ana sürümünün eşleşmesini önerir. Bu nedenle cache `116_GRCh38`
olarak seçildi. `latest` etiketi kullanılmaz.

Kaynak arşiv:

```text
https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
```

Canlı HEAD yanıtında uzunluk **27.644.657.162 bayt** (27,64 GB), Last-Modified
`Thu, 09 Apr 2026 00:45:28 GMT` ve ETag `"66fbffe0a-64efc55395049"` görüldü.
Yanıt başlıkları `reports/generated/local_vep_readiness/preflight.json` içinde
saklandı. HEAD isteği arşivi indirmez; başlıklar arşivin SHA256 doğrulaması yerine
geçmez. Kurulumdan sonra yerel cache dosyaları ayrıca hash ile kilitlenir.

Konteyner kurulduktan sonraki ölçümde C: üzerinde yaklaşık **37,7 GB** boş alan
vardı. Kurulum için belirlenen planlama bütçesi:

```text
2 × arşiv boyutu + 8 GiB çalışma payı = 63.879.248.916 bayt ≈ 63,88 GB
```

Bu değer arşivi saklama, açma ve çalışma payı için ihtiyatlı bir **proje
tahminidir**; sağlayıcının bildirdiği kesin kurulu boyut değildir. Gerçek açılmış
boyut henüz ölçülmedi. Mevcut disk bu bütçeyi karşılamadığı için büyük arşiv
indirmesi başlatılmadı. Docker kurulumu ve VCF hazırlığı tamamlandı.

## Komutlar

Etiketsiz girdiyi hazırlamak veya mevcut dosyaların hash'lerini doğrulamak:

```powershell
.\.venv\Scripts\python.exe scripts\local_vep_pipeline.py prepare
```

Ortam kontrolü; yalnızca resmî arşivin başlıklarını da yenilemek için:

```powershell
.\.venv\Scripts\python.exe scripts\local_vep_pipeline.py check --inspect-remote
```

Hazır olmayan durumda `check` çıkış kodu 2 döndürür ve nedenlerini JSON'a yazar.
Bu koşuda nedenler `release_116_grch38_cache_missing`,
`disk_below_planning_budget`, `resources_not_locked` oldu. Kontrol komutu cache
indirmez; kurulu dosyaların tüm içerik doğrulaması `run` öncesinde yapılır.

Daha geniş disk/sunucu belirlendiğinde `config/local_vep_pipeline.json` içindeki
`cache_dir` o dizine ayarlanabilir. Hazırlık kuralı; girdiler, kapsam ve parça
boyutundan ayrı bir hash üretir. Bu nedenle yalnızca cache/output konumu
değişmesi 170 VCF'yi yeniden oluşturmayı gerektirmez. Kaynak veya parça boyutu
değişirse farklı `prepared_dir` kullanılmalıdır.

Yeterli alana sahip cache dizini oluşturulduktan sonra resmî yükleyici için
doğrulanan seçenekler aşağıdadır. `CACHE_DIZINI` gerçek mutlak yol ile değiştirilir;
hedef dizin Docker bind mount öncesinde mevcut olmalıdır:

```powershell
docker run --rm --mount "type=bind,source=CACHE_DIZINI,target=/cache" ensemblorg/ensembl-vep@sha256:d4f7ab1dad8c350b3ffbda026c1a3086310a6311b53fb90842eaf71c23801def INSTALL.pl --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHE_VERSION 116 --CACHEDIR /cache --CACHEURL https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/ --USE_HTTPS_PROTO
```

Bu kurulum komutu henüz çalıştırılmadı; `INSTALL.pl --help` seçenekleri gerçek
imajda doğrulandı. Modelin mevcut FASTA'sı kullanıldığı için ikinci FASTA
indirilmez. Cache üçüncü taraf veriler içerebilir; önceki
[kaynak/lisans incelemesindeki](context_annotation_pilot.md) sınırlar geçerlidir.

Kurulum tamamlandığında kaynak kilidi:

```powershell
.\.venv\Scripts\python.exe scripts\local_vep_pipeline.py lock-resources
```

Bu adım VEP'in `--show_cache_info` kontrolünü çalıştırır; cache ağacındaki tüm
dosyaların SHA256 envanterini, FASTA/FAI doğrulamasını, imajı ve konfigürasyonu
`resources.lock.json` içine yazar. Bu, yerel dosyaların içerik kilididir;
sağlayıcının imzalı yayın manifesti olarak sunulmaz. Mevcut kilit farklı içerikle
sessizce değiştirilmez. Taşınan veya güncellenen kaynaklar için yeni kilit ve
yeni çıktı dizini kullanılmalıdır.

Önce bir parçalık gerçek çalışma:

```powershell
.\.venv\Scripts\python.exe scripts\local_vep_pipeline.py run --max-chunks 1
```

Kaynaklar hazır olduğunda kalan parçaları işlemek için örnek:

```powershell
.\.venv\Scripts\python.exe scripts\local_vep_pipeline.py run --max-chunks 170
```

Çalışma `--offline`, `--network=none`, `--pull=never` ile yürür. Cache, FASTA ve
girdiler konteynere salt okunur bağlanır. MANE/canonical, sembol, ekzon sırası,
domain ve protein bilgisi istenir; `--check_existing`, patojenite skoru
eklentileri ve klinik veri seçenekleri açılmaz. Ham JSON'da beklenmedik dolu
`colocated_variants` alanı tamamlanma kontrolünde reddedilir.

## Kesintiden devam ve çıktı kontrolü

Her parça kendi deneme dizininde komut, stdout/stderr ve JSONL tutar. Başarılı
çıkış tek başına yeterli değildir: VCF ile tüm kimlikler birebir eşleşmeli;
assembly, gerçek çıktı koordinat/allel alanları, varyant ID'si ve anotasyon
varlığı da doğrulanmalıdır. Atlanan/tekrarlanan kayıt varsa `complete.json`
yazılmaz. Başarısız deneme korunur; sonraki çağrı yeni bir deneme dizininde başlar.

Tamamlanmış parça; girdi, kaynak kilidi, kod/çalışma ayarları ve çıktı hash'leri
aynıysa yeniden hesaplanmaz. Değişiklik varsa işlem durur. Kaynak cache ve FASTA
her `run` başında içerik olarak doğrulanır. Dondurulmuş model ve ClinGen birincil
sonuçlarının hash'leri de işlem öncesi/sonrası korunur. Aynı çıktı dizininde aynı
anda birden çok runner kullanılmamalıdır.

## Doğrulama ve kalan işlem

**10 yeni test geçti.** Gerçek küçük FASTA ile referans kontrolü; train dışlama
ve etiket okumama; parça hesabı; bozulmuş girdi/çıktı/cache reddi; eksik/tekrarlı
varyant kontrolü; yarıda kalan denemeden devam; tamamlanmış parça kullanımı;
çalışma ayarı değişince durma ve cache konumunu değiştirme sınandı.

Gerçek veri hazırlığı ikinci çağrıda mevcut 170 parçanın hash'lerini doğruladı.
Gerçek VEP konteynerinde sürüm, yükleyici yardım seçenekleri ve anotasyon
parametrelerinin ayrıştırılması kontrol edildi. Cache bulunmadığında beklenen
hata kaydedildi. Bu kontroller, 169.265 kaydın anotasyonunun tamamlandığı anlamına
gelmez; cache ile gerçek çalışma ve REST pilotuyla biyolojik çıktı karşılaştırması
henüz yapılmadı.

Hazırlık sonucunun ardından gereken işlem, yeterli disk/sunucu dizininde cache
kurulumu ve ilk gerçek parça çalışmasıdır. Sonrasında yerel transkript geometrisi
kaynağı bağlanmalı ve REST pilotunun splice/protein sütunlarıyla eşleşme
kontrol edilmelidir. Bu runner şimdilik ham VEP JSONL üretir; tam pilot özellik
tablosunu çıkardığı veya modeli geliştirdiği iddia edilmez. Yeni model deneyi
train-only gene-disjoint geliştirme düzeninde yapılmalı; önceki ClinGen ölçümü
yeniden kullanılmamalıdır.
