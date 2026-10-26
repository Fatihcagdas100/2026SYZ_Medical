# ClinGen eRepo kimlik taraması sonucu

Tarih: 26 Ekim 2026 (Europe/Istanbul)

Sonraki aşama tamamlandı: GRCh38 referans ve allel denetimi sonrası 264 kayıt
teknik kabul aldı. Güncel sonuç `docs/clingen_grch38_admission_result.md`
dosyasındadır. Bu belge ilk kimlik taramasının tarihsel kaydıdır.

## Sonuç

Eski 367.196 varyantın VariationID kümesine karşı yapılan taramada **282 aday**
kaldı: 216 Pathogenic ve 66 Benign. Bunlar sadece kimlik düzeyinde eşleşmeyen
kayıtlardır; henüz bağımsız, çakışmasız dış kohort olarak kabul edilmemiştir.

| Tarama aşaması | Sayı |
| --- | ---: |
| İndirilen sınıflandırma kaydı | 13.564 |
| Aktif sınıflandırma kaydı | 13.357 |
| Geri çekilmiş / aktifliği doğrulanmamış kayıt | 207 |
| Aktif kesin Pathogenic/Benign yorum kaydı | 4.826 |
| Bu yorumlarda tek ClinVar kimliği bulunmayan kayıt | 216 |
| Tek ClinVar kimliğine bağlanan kesin ikili varyant | 4.610 |
| Eski veriyle VariationID çakışması nedeniyle dışlanan | 4.321 |
| Aktif yorumları farklı sınıflandırmalar taşıdığı için dışlanan | 7 |
| Eski VariationID kümesiyle eşleşmeyen aday | 282 |

Kimliği olmayan 216 yorum, 'yeni varyant' diye sayılmadı. Bir varyantın farklı
hastalık/panel yorumları varsa bunlardan ilk satır seçilmedi: aktif yorumların
sınıfları farklıysa varyant aday listesinden çıkarıldı. Tutarlı birden fazla
yorumun her biri, hastalık, uzman paneli ve doküman sürümüyle korundu.

## Adayların dağılımı

| Grup | Pathogenic | Benign | Toplam |
| --- | ---: | ---: | ---: |
| Kimliği eşleşmeyen adaylar | 216 | 66 | 282 |
| HGVS başlığına göre tutarlı missense adayları | 142 | 38 | 180 |

Missense sayımı `preferredVarTitle` alanının mevcut HGVS ayrıştırıcısıyla
yorumlanmasına dayanır; bağımsız elle anotasyon denetimi yapılmadı.

## Anlık görüntünün sabitlenmesi

Kaynak, [ClinGen eRepo resmî arayüzünün](https://erepo.clinicalgenome.org/evrepo/)
kullandığı `/evrepo/api/summary/classifications` API'sidir. Toplu TSV indirme
yerine UUID'ye göre sıralanan 28 JSON sayfası yerelde saklandı. İndirme öncesi
ve sonrası toplam sayılar ve kayıt UUID tekilliği doğrulandı. API işlemsel bir
snapshot garantisi vermediğinden, sayım eşitliği indirme sırasında hiçbir
yorumun değişmediğini tek başına kanıtlamaz.

Bu veri, sağlayıcının resmî sürümlü arşivi olarak değil, indirildiği anda
hash'leriyle sabitlenen yerel anlık görüntü olarak tanımlanır. Her yorumun kendi
`docVersion`, onay/yayın tarihi, paneli ve hastalığı ham sayfalarda korunur.
[Kullanım şartlarının](https://www.clinicalgenome.org/docs/terms-of-use/) kopyası
ve hash'i de saklandı.

- Snapshot klasörü: `data/external/clingen_erepo_2026-10-26/`
- Snapshot manifest SHA-256:
  `6d328270433d4f949b4ced5d6e7033e308d11220e4c7b36ea37bebf724a5a5e2`
- Eski veri SHA-256:
  `4d018115be5a7ad710e86fb653630a87a2d9489ad29a59d8cc38b3d9d46341c1`
- Aday yorum tablosu SHA-256:
  `459bae3b90e9da85adcb77dc59bbbca99ad4faabd024cb1fa92ccdd3b2a7c77b`

Ham indirme ve üretilen yerel raporlar sürüm kontrolüne alınmaz. Bu belge ve
çalıştırma betikleri sürüm kontrolünde tutulur.

## Değerlendirme sınırı ve sonraki işlem

282 sayısı aday üst sınırıdır. AlleleID, normalleştirilmiş GRCh38 koordinat-allel
eşdeğerliği, referans doğrulaması ve etiket kökeni henüz incelenmedi. Farklı
ClinVar kimliklerinin aynı biyolojik alleli temsil etmesi mümkündür; ek eleme
sonrası sayı düşebilir.

ClinGen ile ClinVar arasında yorum paylaşımı vardır. Bu yüzden yeni VariationID
bulunması, etiketlerin ClinVar'dan bağımsız üretildiğini göstermez. Kaynak
bağımsızlığı kanıtlanamazsa sonuç uzman-kürasyonlu ek değerlendirme kapsamıyla
tanımlanmalıdır; bağımsız klinik dış doğrulama diye sunulmamalıdır.
[NCBI kaynak açıklaması](https://www.ncbi.nlm.nih.gov/clinvar/intro/).

Sonraki işlem, adayların GRCh38 genom/allel eşlemesini tamamlayıp referanslı kabul
hattından geçirmek ve etiket kökenini belgelemektir. Bu taramada model
yüklenmedi, eşik değiştirilmedi ve dış metrik hesaplanmadı.

## Yeniden üretim

İndirilmiş snapshot'ı tekrar taramak için:

```powershell
.\.venv\Scripts\python.exe scripts\screen_clingen_erepo.py `
  --snapshot-dir data\external\clingen_erepo_2026-10-26
```

Canlı API'yi yeniden indirmek farklı bir anlık görüntü üretir. Eski klasörün
üzerine yazmak yerine yeni bir klasör adı verilmelidir.
