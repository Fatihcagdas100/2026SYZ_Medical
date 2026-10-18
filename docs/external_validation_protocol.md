# Dış doğrulama protokolü

Tarih: 18 Ekim 2026

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

Eşik dış kohortta ayarlanmaz; yeni model eğitilmez; dış metrikler görüldükten
sonra birincil sonucu iyileştirmek için tekrar deneme yapılmaz. Ek analiz gerekirse
ayrı bir keşif analizi olarak işaretlenir ve yeni, bağımsız bir kohortta teyit
edilir.

## Yorum sınırı

Dış doğrulama iyi sonuç verse bile model klinik tanı aracı değildir. ClinVar da
dâhil klinik varyant kaynakları gönüllü bildirimlere ve toplanmış yorumlara
dayanır; model skoru bireysel klinik risk olarak yorumlanamaz. [ClinVar kullanım
notları](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/).
