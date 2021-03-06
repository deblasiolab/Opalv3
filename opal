#!/bin/bash

#Script to identify java VM. Also modifies classpath on cygwin.
#Also used to instantiate headless version of Mesquite.
#Written by Travis Wheeler (www.traviswheeler.com)

#If you find bugs or modify the script to improve functionality, 
# please let us know so we can improve the script for everyone

# Thanks to William Cook (http://www.cs.utexas.edu/users/wcook/) for the bit that
# identifies the path containing this script, and uses it to set the classpath. 
# (the exact site of the work that served as inspiration for this code is lost to antiquity)

#To increase memory allocation adjust the following line:  e.g. to 
#mem="500M" 
mem="1G"


#==============================================
#  Shouldn't need to alter anything below here
#==============================================

# This is intended to perform the job performed by GNU 'readlink -f'
# including following symlinks). The 'readlink' on BSD (and thus OS X)
# doesn't use -f for the same purpose. Credit to G. Nix
# http://stackoverflow.com/a/18442790
realpath_ours() {
  OURPWD=$PWD
  cd "$(dirname "$1")"
  LINK=$(readlink "$(basename "$1")")
  while [ "$LINK" ]; do
    cd "$(dirname "$LINK")"
    LINK=$(readlink "$(basename "$1")")
  done
  REALPATH="$PWD/$(basename "$1")"
  cd "$OURPWD"
  echo "$REALPATH"
}



echo 1>&2

# figure out where I live, so I can run java w/ my containing dir as classpath  
dir=`dirname "$0"`
os=`uname`
if test ${os:0:6} = "CYGWIN"
then
  if test ${dir:1:8} = "cygdrive"
  then
    dir="${dir:10:1}:/${dir:12}"
  fi
fi
#abspath="$(cd "${0%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"
#path_only=`dirname "$abspath"`
path_only=$(dirname $(realpath_ours "$0"))


#================================================================
do_struct=0
make_struct=0
do_facet=0
#determine if structure prediction is required
for arg 
do
   if test $arg = "--use_struct"
   then
       do_struct=1
       make_struct=1
   else
       if test $arg = "--facet" ||  test $arg = "--advising_configuration_file"  ||  test $arg = "--realignment_configuration_file"   ||  test $arg = "--configuration_file"
       then
           do_facet=1
           make_struct=1
           if test $arg = "--advising_configuration_file"  ||  test $arg = "--realignment_configuration_file"  ||  test $arg = "--configuration_file"
           then 
              clean_args[${#clean_args[*]}]=$arg
           fi
       else
           clean_args[${#clean_args[*]}]=$arg
       fi
   fi        
done


#if structure prediction is required, figure out the name of the fasta file
i=0
if test ${make_struct} = 1
then
    for arg in "${clean_args[@]}"
    do
       i=$[$i+1]
       if test $arg = "--in"
       then
           seqfile=${clean_args[$i]}
       fi
       if test $arg = "--in2"
       then
           seqfile2=${clean_args[$i]}
       fi             
    done

    #either the argument "--in" was used, and the sequence file was identified above,
    #or assume it's the final argument
    if [ -z ${seqfile} ]
    then
        seqfile=${!#} 
    fi
    
    if [ ${seqfile} ]
    then    
    
        # create random string, due to http://tldp.org/LDP/abs/html/contributed-scripts.html#PW
        ALPH="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
        LENGTH="16"
        while [ "${n:=1}" -le "$LENGTH" ]
        do
            RNDSTR="$RNDSTR${ALPH:$(($RANDOM%${#ALPH})):1}"
            let n+=1
        done
        
        structfile="${seqfile}.${RNDSTR}.ss"
    
        cmd="time $path_only/predict_structure.pl  --in ${seqfile} --out ${structfile} --tmpdir /tmp"
        echo "Before computing alignment, psipred is used to predict structure, with the following command:" 1>&2 
        echo "$cmd" 1>&2        
        $cmd
        
    fi
    
    if [ ${seqfile2} ]
    then    
    
        # create random string, due to http://tldp.org/LDP/abs/html/contributed-scripts.html#PW
        ALPH="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
        LENGTH="16"
        while [ "${n:=1}" -le "$LENGTH" ]
        do
            RNDSTR="$RNDSTR${ALPH:$(($RANDOM%${#ALPH})):1}"
            let n+=1
        done
        
        structfile2="${seqfile2}.${RNDSTR}.ss"
    
        cmd="time $path_only/predict_structure.pl  --in ${seqfile2} --out ${structfile2} --tmpdir /tmp"
        echo "Before computing alignment, psipred is used to predict structure, with the following command:" 1>&2 
        echo "$cmd" 1>&2        
        $cmd
       
    fi        
fi



#check for --mem flag
i=0
for arg in "${clean_args[@]}"
do
   if test $arg = "--mem"
   then
       mem=${clean_args[$i+1]}
       unset clean_args[$i+1]
       unset clean_args[$i]
   fi
   i=$[$i+1]
done



#===================================================

#figure out where java lives 
if [ $OPAL_JAVA_HOME ]
then
  java="$OPAL_JAVA_HOME/bin/java"
elif [ $JAVA_HOME ]
then
  java="$JAVA_HOME/bin/java"
else
  tmp=`java -version 2>&1`
  if echo "$tmp" | grep -q "command not found"  # no "java", so try "jre"
  then
    tmp=`jre -version 2>&1`
    if echo "$tmp" | grep -q "command not found"
    then
       echo "Can't find java. Try setting either the JAVA_HOME environment variable"
       exit
    else
       java="jre"
    fi
  else
   java="java" 
  fi
fi


#================================================================

#Run Opal.jar with appropriate arguments


cmd="$java -server -Xmx$mem -jar $path_only/bin/Opal.jar"
cmd="$java -server -Xmx$mem -cp $path_only/bin opal.Opal"

if [ ${structfile} ] 
then
    if test ${do_struct} = 1
    then
        cmd="${cmd} --structure_file ${structfile}"
        echo "Now Opal will perform alignment, based on predicted structures" 1>&2
    else
        if test  ${do_facet} = 1
        then
        	cmd="${cmd} --facet_structure ${structfile}"
        	echo "Now Opal will perform alignment and output a Facet score" 1>&2
        fi
    fi
fi 
if [ ${structfile2} ]
then
    cmd="${cmd} --structure_file2 ${structfile2}"
fi   
cmd="${cmd} ${clean_args[@]}"
echo "Running the following java command for Opal:" 1>&2 
echo "$cmd" 1>&2 
$cmd

rm -f ${structfile} ${structfile2}


#When a script ends with an exit that has no parameter, the exit status of the script is the exit status of the last command executed in the script (previous to the exit).
