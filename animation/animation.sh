echo -e "Reissner beam animation"
echo "computing the reissner beam position"
#./reissner --time-resample 256 -s 0.05 --radius 0.1 -o result.tmp --verbose
echo "done"
echo "creating plots... "
#python plot.py result.tmp.group.data result.tmp.energy.data
echo "done\n"
echo "concatenating images..."
for i in $(seq 0 255); do
	plot=`echo img/imgReissner"${i}".png`;
	energy=`echo img/imgEnergy"${i}".png`;
	momenta=`echo img/imgMomenta"${i}".png`;
	convert ${energy} ${momenta} -append tmp.png;
	convert ${plot} tmp.png +append animation/frame${i}.png;
done;
rm tmp.png
echo "done"
echo "creating movie..."
#avconv -i "animation/frame%d.png" -r 25 movie.mov
convert -delay 8 -loop 0 $(for i in $(seq 0 255); do echo animation/frame${i}.png; done) tmp.gif
gifsicle -O3 tmp.gif -o anim.gif --colors 256
rm tmp.gif
echo "done"
