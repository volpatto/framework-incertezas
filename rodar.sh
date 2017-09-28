comando="cd $1"; eval $comando
comando="rm solution*"; eval $comando
comando="rm RC*"; eval $comando
comando="../main.exe | tee tela.txt"; eval $comando
