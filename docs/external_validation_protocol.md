# Dış doğrulama protokolü

Tarih: 18 Ekim 2026

Teknik hat güncellemesi: 26 Ekim 2026

## Amaç

Bu protokol, seçilmiş `HGVS + missense` lojistik modelinin daha önce hiç
kullanılmamış bir kohortta tek birincil değerlendirmesini tanımlar. Dış kohort
geldikten sonra model seçimi, hiperparametre, özellik profili ve karar eşiği
değiştirilmeyecektir.

## Dondurulan aday

`scripts/freeze_external_candidate.py`, yalnızca gene-disjoint geliştirme
ayrımındaki train ve validation satırlarıyla modeli yeniden eğitir. Kilitli iç
test satırları ve dış kohort etiketleri eğitim çerçevesine alınmaz.

- Profil: `hgvs_missense`
- Karar eşiği: `0.50`
- Seed: `2027`
- Artifact ve bütünlük hash'i: `reports/generated/external_candidate/manifest.json`

Bu artifact dış kohort için üretildiğinde henüz bir dış metrik hesaplanmamıştır.
Kohort etiketleri görüldükten sonra yalnızca önceden tanımlı ölçüm yapılır.

## Kohort kabul koşulları

Bir kohort aşağıdaki koşulların tamamını sağlamadan değerlendirmeye alınmaz:

1. Kaynak, sürüm/tarih, indirme adresi, lisans veya kullanım koşulu ve dosya
   SHA-256 değeri kayıtlı olmalıdır.
2. Klinik etiketlerin üretim biçimi ve hedef eşlemesi açık olmalıdır. Bu çalışma
   için yalnızca açıkça `Pathogenic` veya `Benign` olan kayıtlar birincil ikili
   değerlendirmeye girer; belirsiz, `Likely` ve birleşik etiketler ana metrikten
   dışlanır.
3. GRCh38 koordinatı ile referans/alternatif allel veya eşdeğer güvenilir varyant
   tanımlayıcısı bulunmalıdır. HGVS protein gösterimi yoksa modelin eksik-değer
   davranışı raporlanır; eksikliği gizlice doldurulmaz.
4. Geliştirme verisiyle çakışma denetlenir: `VariationID`, `AlleleID` ve
   normalleştirilmiş GRCh38 koordinat-allel anahtarı (`kromozom`, `pozisyon`,
   `referans`, `alternatif`). Çakışan kayıtlar dış değerlendirmeden çıkarılır.
5. Gen çakışması ayrıca raporlanır. Yeterli örnek varsa hiç görülmemiş genler
   için ayrı, ikincil bir gene-disjoint dış alt analiz yapılır.
6. Her iki sınıf da bulunmalı; sınıf sayıları ve prevalans raporlanmalıdır.

### Referans ve normalizasyon

Kabul betiği `--reference-config data/reference/grch38.json` ister. Bu dosya,
kullanılacak **gerçek GRCh38** referansını ve FAI indeksini hash'leriyle sabitler.
Yollar config dosyasının konumuna göre çözülür:

```json
{
  "assembly": "GRCh38",
  "fasta": "GRCh38.fa",
  "fasta_sha256": "DOSYANIN_64_KARAKTER_SHA256_DEGERI",
  "fai_sha256": "FAI_DOSYASININ_64_KARAKTER_SHA256_DEGERI"
}
```

FASTA sıkıştırılmamış, FAI indeksi bu FASTA ile uyumlu olmalıdır. Hash eşleşmesi,
seçilen dosyanın değişmediğini doğrular; assembly'nin kaynağı ve sürümü indirme
kaydından ayrıca doğrulanır. Config'e `GRCh38` yazmak tek başına dosyanın GRCh38
olduğunu kanıtlamaz. Referans 26 Ekim 2026'da UCSC initial hg38 kaynağından
indirilip sağlayıcı MD5 kontrolünden geçirildi. Hazırlama betiği
`scripts/prepare_grch38_reference.py` dosyasındadır. `ChromosomeAccession`
alanını denetlemek için üretilen
`data/reference/ucsc_hg38_initial/grch38_with_aliases.json` kullanılmalıdır.

Her iki kaynaktaki `Assembly` sütunu gereklidir. Kabul algoritması A/C/G/T içeren
iki allelli varyantlarda referans alleli denetler; ortak son/ön ekleri kaldırır ve
indelleri tekrarlı bölgelerde referans bazlarıyla sola hizalar. `chr1`/`1` ve
`chrM`/`MT` aynı kromozom anahtarına çevrilir. Model girdisinde kromozom adı da
uyumlu biçime getirilir. Pozisyon ve allellerden türetilen normalleştirilmiş anahtar
yalnızca çakışma içindir; mevcut modelin allel-uzunluk tarifi değiştirilmez.

Sembolik, çok-allelli, belirsiz bazlı ve eksik alleller dışlanır; bu uygulama
liftover veya karmaşık haplotip eşdeğerliği çözümü yapmaz. Geliştirme verisindeki
normalleştirilemeyen kayıtların kimlikleri yine engel listesinde tutulur ve sayısı
raporlanır. Gen yeniliği tüm eski veriyle karşılaştırılır; bu, iç testte görülmüş
genleri de muhafazakâr biçimde 'görülmüş' kabul eder. Ham test etiketleri üzerinden
metrik hesaplanmaz.

FAI yerleşimi için [HTSlib faidx belgesi](https://www.htslib.org/doc/faidx.html),
normalizasyonun yerleşik araç karşılığı için [BCFtools norm
belgesi](https://samtools.github.io/bcftools/bcftools#norm) kullanılabilir. Bu
repository'nin küçük varyant uygulaması yapay örneklerle sınanmıştır; BCFtools
ile tam uyumluluk iddiası taşımaz.

ClinVar'ın düzenli indirmeleri haftalık güncellenir ve aylık sürümler arşivlenir.
ClinVar'ın daha sonraki bir arşiv sürümünden yalnızca yeni/çakışmasız kayıtlar
oluşturmak **zamansal doğrulama** sağlar; kaynak bağımsızlığı sağlamaz. Bu ayrım
raporda açıkça yazılır. Kaynak bağımsız kohort tercih edilir. [ClinVar indirme
politikası](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/), [ClinVar
sürüm döngüsü](https://www.ncbi.nlm.nih.gov/clinvar/docs/release_cycle/).

## Birincil değerlendirme planı

Kabul edilmiş kohort üzerinde, dondurulmuş artifact ile tek seferde şunlar
raporlanır:

- ROC-AUC ve average precision;
- sabit `0.50` eşiğinde duyarlılık, özgüllük, precision, F1 ve karışıklık matrisi;
- %95 bootstrap güven aralıkları (en az 1.000 tekrar, sabit seed);
- tüm kohort ve missense alt grubu;
- varsa hiç görülmemiş genler alt grubu;
- dış kohortta HGVS/protein anotasyonu eksikliği dağılımı.

Güven aralıkları, sınıf sayılarını koruyan satır bootstrap'ı ile ROC-AUC, AP,
precision, recall, özgüllük ve F1 için hesaplanır. Gen içi bağımlılık bu
aralıklarda modellenmez; yöntem ve bu sınırlama raporda yazılır. Boş veya tek
sınıflı alt gruplar `insufficient_classes` olarak sayılarıyla raporlanır.

### Tek değerlendirme kaydı

Hash ve artifact kontrolleri geçtikten sonra, dış etiketler okunmadan önce
`EVALUATION_STARTED.json` oluşturulur. Aday modelin manifest klasöründeki
`evaluation_registry/<ham-kohort-sha256>.json` kaydı atomik `x` kipinde açılır;
aynı kohortu farklı çıktı klasörüne yönlendirmek de tekrar değerlendirme sağlamaz.
Hata veya yarıda kesilme kaydı silmez. Bu, yerel çalışma disiplini kilididir;
dosyaları elle silmeye veya klasörleri kopyalamaya karşı erişim kontrolü değildir.

Dondurulmuş model üreticisi de dolu artifact/rapor klasörünün üzerine yazmayı
reddeder. Model yeniden eğitilmeden eski artifact korunur.

Eşik dış kohortta ayarlanmaz; yeni model eğitilmez; dış metrikler görüldükten
sonra birincil sonucu iyileştirmek için tekrar deneme yapılmaz. Ek analiz gerekirse
ayrı bir keşif analizi olarak işaretlenir ve yeni, bağımsız bir kohortta teyit
edilir.

## Yorum sınırı

Dış doğrulama iyi sonuç verse bile model klinik tanı aracı değildir. ClinVar da
dâhil klinik varyant kaynakları gönüllü bildirimlere ve toplanmış yorumlara
dayanır; model skoru bireysel klinik risk olarak yorumlanamaz. [ClinVar kullanım
notları](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/).
