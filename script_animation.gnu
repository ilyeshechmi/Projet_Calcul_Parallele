### Animation à partir des fichiers results_t* (x y)

reset

# ---------- Réglages généraux ----------
set term gif animate delay 5 optimize    # GIF animé, delay=5 (modifiable)
set output "animation_results.gif"       # Nom du fichier de sortie

set xlabel "x"
set ylabel "u(x)"
set grid

# Optionnel : fixer les bornes des axes
 set xrange [0:10]
set yrange [0:1.1]

# Récupérer la liste des fichiers results_t* dans l'ordre
files = system("ls resultats/results_t*")

# ---------- Boucle sur tous les fichiers ----------
do for [i=1:words(files)] {
    f = word(files, i)
    print sprintf("Traitement de %s", f)

    # Tracé simple x:y
    plot f using 1:2 with lines title sprintf("File: %s", f)
}

unset output
