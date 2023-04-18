# LMAO
## Lokálna metóda aerodynamického odporu.
## Autor: Radovan Lascsák
## Bakalárska práca, 2023, Univerzita Karlova.
-------------------------------------------

### Abstrakt
Odstraňovanie stredne veľkého kozmického odpadu pomocou lokálne vytvorených oblakov prachových častíc. Na odpad pri prelete oblakom vo vysokej relatívnej rýchlosti pôsobí odporová sila, ktorá ho spomalí a skráti jeho dobu života na orbite. Cieľom je simulácia rozpínania oblakov vo vákuu, preletu odpadu skrz oblak, určenie jeho spomalenia a doby do zhorenia v atmosfére.



### Popis programov:

> **ROS.py** (Rozpínanie Oblaku, Simulácia) - triedy a funkcie na simulovanie rozpínania oblaku vo vákuu a vykresľovanie grafov.
  - **Trieda Oblak:** simulovanie a animovanie relatívnej hustoty oblaku v čase. Jej fitovanie Gausiánom. Fitovanie parametrov Gaussiánu v čase hyperbolou, resp. priamkou.
    - Funkcia colormesh: vytvorenie 3D farebný grafu parametru A Gausiánu po zadanom počte frames (čase) simulácie. Potrebné dáta vygenerované programom oblaky.py.
    - Funkcia colormesh_anim: vytvorenie animácie grafu parametru A v čase. Potrebné dáta vygenerované programom oblaky.py.
  - **Trieda Oblak_BD:** simulovanie relatívnej hustoty oblaku v čase. Ukladanie x_H, t_H. Fitovanie vývoja x_H v čase.
    - Funkcia colormesh_BD: vytváranie 3D farebného grafu stredovej hustoty po zadanom čase t. Vytváranie rôznych grafov zobrazujúcich x_H a t_H. Potrebné dáta vygenerované programom oblaky_BD.py.
    - Funkcia colormesh_BD_anim: vytváranie animácie grafu stredovej hustoty v čase. Potrebné dáta vygenerované programom oblaky_BD.py.
    
> **oblak.py** - vytvorenie jedného oblaku triedy Oblak, alebo Oblak_BD. Práca s ním.

> **oblaky.py** - hromadné spúšťanie simulácií oblakov triedy Oblak multiprocesorovo na zadanom priestore parametrov r a mC. Vytváranie vlastnej priečinkovej štruktúry. Dáta potrebné pre funkcie colormesh a colormesh_anim v ROS-e.

> **oblaky_BD.py** - hromadné spúšťanie simulácií oblakov triedy Oblak_BD multiprocesorovo na zadanom priestore parametrov r a mC. Vytváranie vlastnej priečinkovej štruktúry. Dáta potrebné pre funkcie colormesh_BD a colormesh_BD_anim v ROS-e.

> **oblak_x_H.py** - vykresľovanie závislosti x_H na čase.

> **spracovanie_dat.py** - vykresľovanie grafov pomocou funkcií colormesh, colormesh_anim, colormesh_BD, colormesh_BD_anim v ROS-e.



### Odporúčaný workflow:

Skúmanie správania konkrétneho oblaku (možnosť navolenia r, mC, T, možnosť vykreslenia videa, ukladania, fitovania):
  1. oblak.py - animácia relatívnej hustoty v čase, fit Gausiánom, vykreslenie parametrov Gaussiánu a ich fit v čase
  2. oblak_x_H.py  - závislosť x_H na čase
    
Skúmanie správania veľa oblakov:
  1. generácia dát programom oblaky_BD.py (navolenie T a priestoru r,mC)
  2. spracovanie_dat.py (potrebné zadať rozmery priestoru r,mC a teplotu, navolenie skúmaného času t) - colormesh_BD na grafy, colormesh_BD_anim na animáciu
    
 
### Plánované additions:
 > POS.py (Prelet Oblakom, Simulácia) - zmena rýchlosti odpadu pri prelete profilom hustoty oblaku zisteného pomocou ROS-u.
