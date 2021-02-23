# Svhip - Global ToDo

Sortiert nach interner Prioritaet! Erledigte Punkte werden mit (+) markiert. 

## Refactor Readme.md
Die aktuelle Readme beschreibt die Core-Funktionalitaet, aber geht nicht weiter auf Feinheiten der Funktionsparameter, einlesen multipler Alignments etc. ein. Folgende Aspekte sollten noch abgedeckt werden:
Kommandos fuer das einlesen mehrerer Alignments, unterbinden der automatischen Negativset-Generierung, Verweis auf ein konkretes Beispiel (dem Repo hinzufuegen!).

## Provide a working example (+)
Da bereits Beispielordner existieren, sollte ein ausgewaehltes Example dem Repo hinzugefuegt werden, s. o. 


## Include RNAzWindow.pl Script

Aktuell erfordert Svhip eine funktionale RNAz Installation. Mit eurem Einverstaendnis wuerde ich das rnazWindow script ebenfalls inkludieren und die Bindings im Code direkt loesen anstatt eine Commandline aufzurufen. Dies duerfte die stabilere Loesung sein. 

## Delete unnecessary PyCache
Es wurden am Anfang irrtuemlich einige .pycache Dateien geadded. Diese sind zu entfernen.

## Refactor function parameters
Manche Funktionsdefinitionen sind in ihrem Umfang zu ausladend. **kwargs-Konstruktionen und ggf Subfunktionen sollten genutzt werden um den Code weiter zu entschlacken. (--> lang) 

## Refactor file structure / delete redundant
Trotz der letzten Ueberarbeitung verbleiben aufgrund verschobener Anforderungen einige unnoetig komplizierte Ordnerstrukturen / unintutive Dateisortierung & Funktionsplatzierung. Das betrifft die Entwicklerseite mehr als die Nutzerseite, sollte aber noch mal angegangen werden. (--> lang) 

## Set-basierte Implementierung von Objektsammlungen
Set-Operationen sind im allgemeinen in Python schneller als Listen-/Arrayoperationen. Daher sollten fuer eine Runtime-optimierte Herangehensweise Sets & Hashmaps bevorzugt werden wo immer die Reihenfolge keine Rolle spielt (Sequenzvergleiche, Window-Verwaltung etc...). Keine funktionale Prioritaet sondern reine Optimierung, daher am Ende der Liste. 
