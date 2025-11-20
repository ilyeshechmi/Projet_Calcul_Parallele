reset

# ---------- Réglages généraux ----------
set term gif animate delay 5 optimize
set output "animation.gif"

set xlabel "x"
set ylabel "u(x)"
set grid

# Optionnel : fixer les bornes des axes
set xrange [0:10]
set yrange [-1:2]

# Récupérer la liste des fichiers results_t* triés par le temps (numérique)
# - printf '%s\n' ... pour avoir un fichier par ligne
# - sort -t_ -k2.2n : sépare au "_" et trie numériquement à partir du 2ème caractère de "t0.000000.dat"
files = system("printf '%s\n' resultats/results_t* | sort -t_ -k2.2n")

nfiles = words(files)

# ---------- Boucle sur tous les fichiers ----------
do for [i=1:nfiles] {
    f = word(files, i)
    print sprintf("Traitement de %s", f)

    # Tracé simple x:y
    plot f using 1:2 with lines title sprintf("File: %s", f)
}

unset output
