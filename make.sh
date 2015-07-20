javac facet/*.java opal/*.java opal/*/*.java opal/*/*/*.java
exit;

cd facet
javac *.java
cd ..

cd opal
sh make.sh
cd ..

