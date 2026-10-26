# ClinGen GRCh38 eşleme ve teknik kabul sonucu

Tarih: 26 Ekim 2026

Takip ölçümü tamamlandı: `docs/clingen_supplementary_evaluation_result.md`.
Bu belge model puanlanmadan önce yapılan teknik kabulün tarihsel kaydıdır.

## Sonuç

İlk kimlik taramasındaki 282 adaydan **264 varyant** referans ve teknik çakışma
denetiminden geçti: 198 Pathogenic ve 66 Benign. Bu, uzman-kürasyonlu ek
değerlendirme kohortudur (`expert_curated_supplementary`). Kaynak bağımsızlığı
doğrulanmamıştır; hiçbir model skoru hesaplanmadı.

| Aşama | Sayı |
| --- | ---: |
| İlk eRepo kimlik adayı | 282 |
| Eski yerel ClinVar GRCh38/germline kaydına tekil bağlanan | 156 |
| Registry ve referans ile doğrulanan | 264 |
| Geçerli CAID olmadığı için çözülemeyen | 16 |
| Registry–ClinVar kimlik bağı doğrulanamadığı için dışlanan | 2 |
| 264 kayıtta ek VariationID/AlleleID/normalleştirilmiş koordinat çakışması | 0 |
| Son teknik kabul | 264 |

126 aday eski yerel ClinVar dosyasında bulunamadığı için yalnızca yerel kayıtla
eşleme yeterli değildi. ClinGen Allele Registry'den 266 geçerli CAID için yanıt
alındı; bunların hash ve indirme zamanı kayıtları yerelde korundu. CAID olmayan
16 kayıt belirsiz/karmaşık yapısal gösterimler de içeriyor; bunlar küçük iki
allelli varyant hattına zorla dönüştürülmedi.

## Son anotasyon dağılımı

| Grup | Pathogenic | Benign | Toplam |
| --- | ---: | ---: | ---: |
| Kabul edilmiş kohort | 198 | 66 | 264 |
| Missense | 148 | 39 | 187 |

32 kayıtta protein anotasyonu yoktur. İlgili genlerin tamamı eski veride
görülmüştür: görülmemiş gen kaydı **0**. Bu kohortla yeni genlere genelleme
iddiası kurulamaz. Pathogenic prevalansı `%75` olduğundan AP yorumlanırken bu
taban oranı dikkate alınmalıdır.

İlk taramadaki 180 missense ile son 187 missense sayısının farkı anotasyon
kaynağı değişiminden gelir. Name alanı 154 kayıtta sabitlenmiş eski ClinVar'dan,
109 kayıtta eRepo tercih edilen başlığından, 1 kayıtta Registry standart
başlığından alındı. Bu seçim kuralı skor görülmeden uygulanmıştır; anotasyon
kaynakları satır bazında saklanır.

## Referans ve koordinat denetimi

[UCSC initial hg38 / GCA_000001405.15](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
FASTA'sı indirildi. Sıkıştırılmış dosya MD5 değeri sağlayıcının
[resmî sağlama toplamıyla](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt)
doğrulandı: `1c9dcaddfa41027f17cd8f7a82c7293b`.

FAI indeksi oluşturulup kromozom uzunlukları sağlayıcının `hg38.chrom.sizes`
dosyasıyla karşılaştırıldı. FASTA, FAI ve `chromAlias` SHA-256 değerleri
`data/reference/ucsc_hg38_initial/grch38_with_aliases.json` içinde sabitlendi.

Gerçek ClinVar kontrolü, `Chromosome=9` olup `ChromosomeAccession=NW_...` taşıyan
alternatif dizi kayıtlarını ortaya çıkardı. Bu koordinatlar artık ana kromozom
koordinatı sayılmıyor. Registry dönüşümünde yalnızca doğrulanan GRCh38 ana
kromozom NC_ erişim numaraları kullanılıyor. Registry'nin 0 tabanlı aralıkları,
referanstan gerekli çapa bazı eklenerek VCF gösterimine çevriliyor; ardından REF
eşleşmesi ve indel sola hizalaması denetleniyor.

Eski 367.196 kaydın 5.713'ünde koordinat normalizasyonu uygulanamadı: 5.711
desteklenmeyen/eksik allel ve 2 ana kromozoma eşlenmeyen accession. Bu kayıtların
VariationID ve AlleleID değerleri engel kümesinde tutuldu. 'Koordinat çakışması
yok' bulgusu normalleştirilebilen küçük varyant kapsamıyla sınırlıdır.

## Etiket kökeni

İkili hedefler eRepo uzman yorumlarından gelir; eski ClinVar birleşik/çelişen
etiketleri hedefe çevrilmedi. Her satırda eRepo yorum UUID'si, hastalık, panel,
sürüm ve varsa eski ClinVar etiketi tutulur. Referans/allel bilgisi Registry'den
alınması etiketlere kaynak bağımsızlığı kazandırmaz. Bu nedenle teknik kabul
manifestinde `source_independence_verified=false` korunmuştur.

## Yerel çıktılar

- Kabul manifesti: `reports/generated/clingen_reference_admission/manifest.json`
- Kohort: `data/external/clingen_reference_checked.tsv.gz`
- Özellik tablosu: `data/external/clingen_missense_features.tsv.gz`
- Özellik manifesti: `reports/generated/clingen_external_features/manifest.json`
- Çözülemeyen kayıt denetimi:
  `reports/generated/clingen_allele_resolution_v2/resolution_audit.tsv`

Kabul kohortu SHA-256:
`fea76fa8c6107b8d20e95b50d22a68aaac45cb84b16686018b31e04d7cd897e6`

Model eğitilmedi veya puanlanmadı. Sonraki değerlendirmede dondurulmuş aday ve
sabit `0.50` eşik kullanılmalı; sonuçlar bu kohortun ek değerlendirme kapsamıyla
raporlanmalıdır. Eski iç testte yeni ölçüm yapılmadı.
