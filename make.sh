javac -source 1.5 -target 1.5 facet/*.java opal/*.java opal/*/*.java opal/*/*/*.java
exit;

cd facet
javac -source 1.5 *.java
cd ..

cd opal
sh make.sh
cd ..

