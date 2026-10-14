for f in `find . -type f -name "equi.gsd"`
do
  if [ ! -f $(dirname $f)/polymerize.gsd ]; then
    echo $f
#    echo $(wc -c <"$f")
    rm $f
#    rm $(dirname $f)/polymerize.txt
  fi
  # mv $f $(dirname $f)/deform.gsd
done
