# Impostazioni generali
set terminal gif animate delay 5 loop 0  # Genera una GIF animata con un ritardo di 10ms tra i frame
set output 'animation.gif'                # Nome del file GIF di output

set xlabel 'x'                             # Etichetta asse X
set ylabel 'Function values'               # Etichetta asse Y
set title 'Time Evolution of Function'     # Titolo del grafico

# Trova il numero di colonne nel file
stats 'Animation_points.txt' nooutput
num_cols = STATS_columns

# Loop su tutte le colonne eccetto la prima (da 2 a num_cols)
do for [i=2:num_cols] {
	set xrange [-3:3]
	set yrange [-0.5:3.5]
    plot 'Animation_points.txt' using 1:i with lines title sprintf('Step %d', i)
}
unset output
print "GIF creata con successo: animation.gif"

