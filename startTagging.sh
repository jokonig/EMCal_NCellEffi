if [ $1 == "nomB" ]; then
  root -q -b 'MakeNCellEffiPureGammas.cxx+(1, "13TeVNomB", 0, 0)' &
  root -q -b 'MakeNCellEffiPureGammas.cxx+(0, "13TeVNomB", 0, 0)'
  wait

elif [ $1 == "light" ]; then
  root -q -b 'MakeNCellEffiPureGammas.cxx+(1, "13TeVLowB", true, 20)' &
  root -q -b 'MakeNCellEffiPureGammas.cxx+(0, "13TeVLowB", true, 20)'
  wait
elif [ $1 == "lowB" ]; then
  root -q -b 'MakeNCellEffiPureGammas.cxx+(1, "13TeVLowB", 0, 0)' &
  root -q -b 'MakeNCellEffiPureGammas.cxx+(0, "13TeVLowB", 0, 0)'
  wait
fi
