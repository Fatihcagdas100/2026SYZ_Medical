# Splice ve protein bağlamı: anotasyon uygulanabilirlik pilotu

Bu adım, önceki hata analizinden sonra seçilen özellik kaynaklarını küçük bir
örneklemde sınar. Amaç ekzon sınırı, transkript ve protein bölgesi bilgisinin
GRCh38 varyantlarına doğrulanabilir biçimde eklenip eklenemediğini görmek.
Yeni model eğitimi veya başarı değerlendirmesi değildir.

## Kaynak seçimi ve kullanım koşulları

26 Ekim 2026 tarihinde incelenen birincil kaynaklar:

| Kaynak | Sağladığı bilgi | Bu adımda karar |
| --- | --- | --- |
| [Ensembl VEP REST](https://rest.ensembl.org/documentation/info/vep_region_get) | Transkript başına sonuç terimleri, MANE/canonical işaretleri, ekzon sırası, protein konumu, örtüşen domain anotasyonları | Küçük pilot için kullanıldı; beklenen Ensembl veri sürümü 116 olarak kilitlendi. |
| [Ensembl transcript lookup](https://rest.ensembl.org/documentation/info/lookup) | Sürümlü transkript, ekzon koordinatları, strand ve protein uzunluğu | Gerçek splice birleşimlerine uzaklık hesaplamak için kullanıldı. |
| [Ensembl kullanım koşulları](https://jun2026.archive.ensembl.org/info/about/legal/disclaimer.html) | Projenin ürettiği veri için kullanım kısıtı uygulanmadığı; proje kodunun Apache 2.0 olduğu; üçüncü taraf kısıtlarının ayrıca geçerli olabileceği açıklanıyor | Özellikle Pfam/InterPro gibi dış kaynakların koşulları, toplu dağıtım veya ürünleştirme aşamasında ayrıca kaydedilmeli. Tüm VEP çıktısının tek bir lisansa sahip olduğu varsayılmadı. |
| [Broad SpliceAI Lookup](https://spliceailookup.broadinstitute.org/) | Dizi tabanlı splice tahmini; kullanılan transkript/model ayarları sonuçları etkiliyor | Site kullanıcı başına dakikada yalnızca birkaç sorgu desteklediğini, toplu işler için yerel kurulum önerdiğini belirtiyor. Toplu API taraması yapılmadı. Site SpliceAI için akademik/ticari olmayan kullanım ve GPLv3/CC BY-NC 4.0 koşullarını açıklıyor. |

Ensembl ekzon geometrisi bir SpliceAI skoru değildir. Burada yeni bir splice
tahmin modeli çalıştırılmadı. Başka araçların patojenite skorları, ClinVar
etiketleri, uzman gerekçeleri veya hastalık sınıflandırmaları özellik yapılmadı.
Genin hastalık mekanizması da bu pilotta çözülmüş sayılmaz.

## Örnekleme ve veri sınırı

`config/context_annotation_pilot.json` giriş dosyalarının SHA256 değerlerini,
örnekleme tohumunu, altı tabakayı ve iki örnek/tabaka kuralını kilitler.
Yalnızca **train** bölümünden 12 birincil kromozom SNV'si seçilir: missense,
synonymous, stop-gained/changed, donor +1/+2, acceptor -1/-2 ve deep-intronic.
Tabakalar mevcut HGVS özelliklerinden gelir. Her tabakada `SHA256(seed:VariationID)`
sırası kullanılır; etiket veya model tahmini örneklemeyi etkilemez. Küçük ve
tabakalı örneklemden tüm veri seti için kapsam oranı çıkarılamaz. İndeller,
yapısal varyantlar, alternatif contigler ve diğer varyant grupları sınanmadı.

Daha önce hata analizinde görülen HBB `38682`, ACADVL `569548` ve GFAP `190351`
ayrıca **önceden görülmüş tanısal örnek** olarak anotasyonlanır. Bunlar train
özellik tablosuna katılmaz. Sonuçlarının incelenmesi yeni bağımsız doğrulama
değildir; klinik etiketleri ve önceki olasılıkları bu betik tarafından okunmaz.

## Üretilen sütunların anlamı

| Sütun | Anlamı |
| --- | --- |
| `biotype` | Seçilen transkriptin protein-coding vb. türü. |
| `consequence_terms` | Seçilen transkript için VEP sonuç terimleri; aynı varyant hem synonymous hem splice-region olabilir. |
| `exon_rank`, `exon_total` | Varyantın bulunduğu ekzonun transkript içindeki sırası ve toplam ekzon sayısı. İntronik varyantlarda boş kalır. |
| `location` | Seçilen transkripte göre exonic, intronic veya outside_transcript. |
| `donor_distance_bp`, `acceptor_distance_bp` | İlgili gerçek splice birleşimine baz cinsinden geometrik uzaklık; strand dikkate alınır. Ekzonun terminal bazı için 0, ilk intronik baz için 1'dir. |
| `nearest_splice_boundary_bp` | Tanımlı donor/acceptor uzaklıklarının küçüğü. Transkriptin 5'/3' uçları splice birleşimi sayılmaz; bulunmayan birleşim boş kalır. |
| `protein_relative_position` | VEP protein başlangıç konumu / seçilen proteinin uzunluğu. Protein sonucu bulunmadığında boş kalır. |
| `pfam_annotation_count`, `interpro_annotation_count` | Varyantla örtüşen benzersiz Pfam/InterPro anotasyonu sayısı. Sıfır, bu yanıtta anotasyon bulunmadığını gösterir; işlevsel etki yokluğu anlamına gelmez. |

Tek transkript seçimi etikete bakmadan şu sabit sırayla yapılır: kaynak gen
sembolüyle eşleşme, MANE Select, canonical, transkript kimliğinin alfabetik sırası.
Şiddet veya klinik sonuç bu sırayı değiştirmez. Alternatif transkriptler ham
yanıtlarda korunur; `annotation_audit.tsv` farklı sonuç kümelerini gösterir.
Bu kural bir transkript seçme politikasıdır; her hastalık için doğru izoformu
bulduğuna ilişkin iddia değildir.

`VariationID` yalnızca birleştirme anahtarıdır. Gen/transkript kimlikleri ve tam
koordinatlar denetim içindir; aday biyolojik özellik listesine alınmaz.

## Çalıştırma ve yeniden üretim

İlk indirme için:

```powershell
.\.venv\Scripts\python.exe scripts\pilot_context_annotations.py --fetch
```

Hash denetimli yerel yanıtlarla, ağ kullanmadan yeniden üretmek için:

```powershell
.\.venv\Scripts\python.exe scripts\pilot_context_annotations.py
```

Betik seçilen varyantları sabit GRCh38 FASTA ile doğrular. VEP yanıtındaki
assembly, kromozom, konum, REF/ALT ve lookup transkript sürümü/strand eşleşmesi
zorunludur. REST yanıtları URL, UTC zamanı ve SHA256 makbuzlarıyla
`data/external/ensembl_context_pilot_2026-10-26/` altında korunur. İndirme
modunda canlı sürüm başlangıçta ve sonda ayrıca denetlenir; eski önbelleğe yeni
sürümden yanıt eklenmesi reddedilir. Ağsız yeniden üretim yalnızca korunmuş
yanıtların sürümünü doğrular; canlı servisin değişmediği iddiasında bulunmaz.

Bu bir yerel erişim anı görüntüsüdür. REST verisinin resmî, indirilebilir VEP
cache sürümüyle birebir aynı olduğu varsayılmaz. Büyük ölçekli sonraki adımda
yerel VEP/cache/FASTA sürümleri ve tüm dış kaynaklar ayrıca kilitlenmelidir.

`reports/generated/context_annotation_pilot/` çıktıları:

- `annotation_audit.tsv`: 12 eğitim ve 3 önceden görülmüş örneğin seçimi ve anotasyon ayrıntıları.
- `train_candidate_features.tsv`: yalnızca eğitim örnekleri için izin verilen biyolojik sütunlar.
- `seen_diagnostic_features.tsv`: önceden görülmüş üç örneğin ayrı özellik tablosu.
- `annotation_failures.tsv`: yeniden denemelerden sonra anotasyon alınamayan kayıtlar; başarı paydasından gizlenmez. Korunmuş hata makbuzları ağsız modda da yeniden okunur.
- `manifest.json`: giriş/çıktı hash'leri, referans ve servis sürümleri, yanıt makbuzları ve kapsam sayıları.
- `fetch_manifest.json`, `replay_manifest.json`: son canlı sürüm denetimli erişim ve son ağsız yeniden üretim kayıtları ayrı saklanır.

Dondurulmuş model, ClinGen ölçüm manifesti ve kaydedilmiş tahmin dosyasının
hash'leri başlangıçta ve sonda doğrulanır. Model yükleme, tahmin üretme, eşik
arama veya metrik yeniden hesaplama yapılmaz.

## 26 Ekim 2026 pilot sonucu

**15/15 varyant başarıyla anotasyonlandı:** 12 train ve 3 önceden görülmüş tanısal
örnek. BLNK `1220868` ilk indirmede üç denemede de HTTP 500 yanıtı verdi;
önbellekten devam eden sonraki koşuda aynı varyant ve aynı parametrelerle
anotasyon alındı. Örnek değiştirilmedi. Tamamlanan koşunun hata tablosu boştur.

| Kontrol / anotasyon | Train örnekleri | Önceden görülmüş örnekler |
| --- | ---: | ---: |
| Denenen / tamamlanan | 12 / 12 | 3 / 3 |
| Kaynak gen sembolü eşleşen | 12 / 12 | 3 / 3 |
| MANE Select seçilen | 12 / 12 | 3 / 3 |
| Gerçek splice birleşimi uzaklığı hesaplanan | 12 / 12 | 3 / 3 |
| Protein göreli konumu bulunan | 6 / 12 | 3 / 3 |
| Pfam örtüşme anotasyonu bulunan | 2 / 12 | 3 / 3 |
| Yanıtta doğrudan `InterPro` kaydı bulunan | 0 / 12 | 0 / 3 |
| Transkriptler arasında birden fazla sonuç kümesi bulunan | 10 / 12 | 3 / 3 |

Protein konumu bulunmayan altı train örneği introniktir. InterPro sütunu bu
örneklemde bilgi taşımadı; VEP yanıtında bu isimle kayıt bulunmaması InterPro'da
biyolojik anotasyon olmadığı anlamına gelmez. Pfam kapsamı da tüm varyantlara
genellenemez. Alternatif transkriptlerin farklı sonuç vermesi, seçim politikasını
ve ham yanıtların korunmasını gerekli kılar.

Üç tanısal örnekte gözlenen bağlam:

| Varyant | Seçilen transkript | VEP sonucu | Splice / protein bağlamı |
| --- | --- | --- | --- |
| HBB `38682` | ENST00000335295.4 | synonymous + splice-region | Ekzon 1/3; gerçek donor sınırına ekzon içinden 2 baz uzaklık. |
| ACADVL `569548` | ENST00000356839.10 | synonymous + splice-region | Ekzon 12/20; donor sınırındaki son ekzon bazı, uzaklık 0. |
| GFAP `190351` | ENST00000588735.3 | stop-gained | Ekzon 6/9; göreli protein konumu yaklaşık 0.722. Bu bilgi tek başına hastalığın fonksiyon kaybı/kazanımı mekanizmasını veya uzman Benign yorumunu açıklamaz. |

HBB ve ACADVL için önceki kaba `synonymous` gösterimin ötesinde ölçülebilir
bağlam elde edildi. Bunun duyarlılığı artırdığı henüz gösterilmedi; model
yeniden eğitilmedi ve bu örnekler yeni doğrulama sayılmadı. SpliceAI skoru,
deneysel RNA sonucu, NMD veya hastalık mekanizması bu pilotun çıktısı değildir.

REST API sürümü `15.12`, Ensembl yazılım/veri sürümü `116` olarak kaydedildi.
Gerçek önbellek üzerinde ağ çağrıları hata verecek şekilde engellenerek yeniden
üretim çalıştırıldı: dört TSV dosyasının SHA256 değeri canlı erişim koşusuyla
birebir aynı kaldı. Dondurulmuş model ve birincil ClinGen sonuç hash'leri korundu.
Yeni betiğin **9 testi** geçti; örnekleme, iki strand üzerindeki uzaklıklar,
etiket/skor dışlama, koordinat ve transkript sürümü reddi, hata makbuzu,
bağımsız hata tablosu ve sürüm yükselmesinde durma kontrol edildi.

### Bu koşunun bütünlük kaydı

| Dosya | SHA256 |
| --- | --- |
| `fetch_manifest.json` | `bc005a1976201fe6979144d08e5aa8aba91d0d3af6359d3eace4d6a5edc65bcb` |
| `replay_manifest.json` | `81d9a3819bd435091e4a6715666224e4af2bc7617bcd8eaae22706d6ee9ac5a6` |
| `train_candidate_features.tsv` | `699fea68a4912b1ad4563775b550f4f6027382b47ff2f39d8b0ccdc0c3a2890b` |
| `seen_diagnostic_features.tsv` | `4fc57480a2509e073434e00eb39f9410206dfa28a8c22626d95a42dd249b8791` |
| `annotation_audit.tsv` | `01374805596399ed1d6ae81266f61171b5593bf300da04bc0d95f09cbf7db297` |

Sonraki işlem, REST pilotunu büyük veri için yerel ve sürümü sabit bir VEP
işlem hattına taşımaya hazırlanmaktır. Öncelik splice sınır geometrisi ve
transkript/protein bağlamıdır; sürekli sıfır kalan InterPro sütunu bir gelişme
iddiasını desteklemez. Ardından train-only gene-disjoint geliştirme deneyinde
mevcut özellikle yeni bağlamlı gösterim karşılaştırılabilir. Bir sonraki modelin
dış başarısı için mevcut ClinGen kohortundan ayrı, kullanılmamış veri gerekir.
