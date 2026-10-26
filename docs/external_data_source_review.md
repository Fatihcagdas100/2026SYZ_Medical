# Açık dış veri kaynaklarının ön incelemesi

İnceleme tarihi: 26 Ekim 2026

## Karar

İlk teknik uygunluk taraması için **ClinGen Evidence Repository** adayı öne
alınmıştır: uzman klinik sınıfları ile varyant kimlikleri aynı arayüzde sunulur.
Bu bir **kaynak-bağımsızlık onayı değildir**. Sürüm sabitlenmesi, etiket kökeni,
kullanım şartları ve eski veriyle çakışma sayıları incelenmeden bu adaydan dış
performans raporlanmaz.

Ön incelemenin ardından ClinGen eRepo anlık görüntüsü indirilip VariationID
çakışması tarandı. 282 kimliği eşleşmeyen aday kaldı; allel/koordinat ve
referans denetiminin ardından 264 kayıt uzman-kürasyonlu ek kohort olarak teknik
kabul aldı. Kaynak bağımsızlığı doğrulanmadı.
Kabul sonucu `docs/clingen_grch38_admission_result.md` dosyasındadır. Sonraki tek
ölçümde genel ROC-AUC `0.7932`, missense ROC-AUC `0.6642` bulundu; ayrıntılar
`docs/clingen_supplementary_evaluation_result.md` dosyasındadır. Aşağıdaki kaynak
uygunluğu yorumları dokümantasyondan yaptığımız çıkarımlardır.

## Kaynak karşılaştırması

| Kaynak | Etiket / ölçüm | Erişim | Projemiz açısından karar |
| --- | --- | --- | --- |
| ClinGen Evidence Repository | Uzman-kürasyonlu patojenite sınıfları | TSV/CSV indirme seçenekleri ve API dokümanı | Klinik ikili hedef için ilk tarama adayı; ClinVar kimliği ve kaynak ilişkisi incelenmeli |
| BRCA Exchange / ENIGMA | BRCA1/2 uzman sınıfları ve birleştirilmiş kaynak bilgisi | Sürümlü veriler ve kaynak raporlarını ayırabilen API | Gen kapsamı dar; ClinVar bir kaynak olduğundan bağımsızlık varsayılamaz |
| MaveDB / MaveMD | Deneysel fonksiyon skorları ve kalibrasyonlar | Tek score-set CSV/API; arşiv seçeneği | Missense için ayrı fonksiyonel doğrulama adayı; mevcut klinik ikili hatta doğrudan sokulmaz |

### ClinGen Evidence Repository

Resmi arayüz uzman sınıflandırmalarını ve destekleyici kanıt özetlerini sunar;
ClinVar Variation ID alanı ile TSV/CSV indirme seçenekleri bulunur.
[Resmi kaynak](https://erepo.clinicalgenome.org/evrepo/), [API
dokümanı](https://erepo.genome.network/docs/cg-erepo/).

ClinVar, ClinGen uzman panellerinin yorumlarını da arşivler. Dolayısıyla başka
bir web sitesinden indirilmiş olmak etiketin bağımsız olduğu anlamına gelmez.
[NCBI kaynak açıklaması](https://www.ncbi.nlm.nih.gov/clinvar/intro/).

Tarama sırasında panel, hastalık, sınıflandırma tarihi, sınıflandırma sürümü,
kanıt kökeni ve aynı varyant için çelişen yorumlar korunmalı; sonradan salt
`VariationID` üzerinden ilk satır seçilmemelidir. `Likely`/VUS kayıtları mevcut
kesin ikili hedefe çevrilmez. [Kullanım şartları](https://www.clinicalgenome.org/docs/terms-of-use/)
indirme anında kaydedilmelidir.

### BRCA Exchange / ENIGMA

Resmi açıklamaya göre platform ClinVar, LOVD, gnomAD ve başka kaynakları
birleştirir. Kaynak raporları API üzerinden ayrıştırılabilir. Bu nedenle yalnızca
ENIGMA kökeni belirlenebilen klinik yorumlar aday olabilir; birleşik tabloyu
bağımsız etiket kaynağı gibi kullanmayacağız. İki genle sınırlı sonuç bütün
genlere genelleme iddiası taşımaz.
[Kaynak ve birleştirme yöntemi](https://www.brcaexchange.org/trackhubs/description.html),
[API kaynak raporları](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/content/api_docs/api_overview.md).

Sürüm ve veri kullanım şartlarının yanı sıra, ENIGMA yorumunun eski ClinVar
eğitim kaydında bulunup bulunmadığı kayıt düzeyinde araştırılmalıdır.

### MaveDB / MaveMD

Ham deney skorları, ClinVar klinik kontrol anotasyonlarından ayrılabilir.
MaveMD kalibrasyonlarında ClinVar kontrol varyantları da kullanılır; bu kontrol
etiketlerini yeniden 'bağımsız dış hedef' diye kullanmayacağız.
[İndirme biçimleri](https://www.mavedb.org/docs/mavedb/finding-data/downloading.html),
[ClinVar entegrasyonu](https://www.mavedb.org/docs/mavedb/mavemd/index.html).

Bir deneyin normal/anormal fonksiyon ölçümü, projemizin Pathogenic/Benign klinik
etiketiyle aynı hedef değildir. Böyle bir çalışma için ayrı fonksiyonel
değerlendirme protokolü gerekir: deneyde skorun yönü, yayımlanmış eşikler,
GRCh38 eşleşmesi ve eski varyantların dışlanması sonuç görülmeden sabitlenir.
Kullanım lisansı score-set bazında kontrol edilir; tüm veri için tek bir lisans
varsayılmaz. [Lisans açıklaması](https://www.mavedb.org/docs/mavedb/faq.html#licensing).

## İlk kohortu seçme ölçütü

Kaynak önceliği model skoruna göre değiştirilmeyecektir. Seçim için şunları
ölçeceğiz: kaynak/sürüm/lisansın sabitlenebilir olması, hedefin klinik anlamı,
GRCh38 ve HGVS kapsamı, eski veriyle çakışma sonrası kalan iki sınıfın sayıları,
missense sayısı ve etiketlerin ClinVar'dan aktarılmış olup olmaması.

Şimdiki dış hat teknik kabulü otomatikleştirir; `--cohort-kind independent`
argümanı kaynak kökeni hakkında kullanıcının beyanıdır ve tek başına bağımsızlığı
kanıtlamaz. Bu ön inceleme sonucunda henüz kaynak-bağımsız klinik dış kohort
seçilmemiştir. Yalnızca zamansal ya da uzman-kürasyonlu ek değerlendirme
sağlanabiliyorsa sonuç o kapsamla raporlanacaktır.
