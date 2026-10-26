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
- `manifest.json`: giriş/çıktı hash'leri, referans ve servis sürümleri, yanıt makbuzları ve kapsam sayıları.

Dondurulmuş model, ClinGen ölçüm manifesti ve kaydedilmiş tahmin dosyasının
hash'leri başlangıçta ve sonda doğrulanır. Model yükleme, tahmin üretme, eşik
arama veya metrik yeniden hesaplama yapılmaz.
