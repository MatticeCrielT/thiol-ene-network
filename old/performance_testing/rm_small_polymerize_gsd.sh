for f in `find . -type f -name "polymerize.gsd"`
do
  if [ $(wc -c <"$f") -le 6000 ]; then
    echo $f
    echo $(wc -c <"$f")
    rm $f
    rm $(dirname $f)/polymerize.txt
  fi
  # mv $f $(dirname $f)/deform.gsd
done
