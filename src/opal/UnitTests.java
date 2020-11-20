package opal;

import com.traviswheeler.libs.*;

import opal.IO.*;
import opal.makers.*;
import opal.tree.Tree;
import opal.exceptions.GenericOpalException;
import facet.*;
import opal.realignment.realignmentDriver;

class UnitTests{
  public static void main(String[] argv){
    test_argHandler_defaults_dna();
    test_argHandler_defaults_protein();
    test_argHandler_defaults_protein_2files();

    System.out.println("All tests passed!");
  }
  
  static void test_argHandler_defaults_dna(){
    String[] argv = {"dna_2seq.fa"};
    ArgumentHandler argHandler = new ArgumentHandler(argv);
    assert(argHandler.getFileA() != null);
    assert(argHandler.getFileB() == null);
    assert(argHandler.getStructFileA() == null);
    assert(argHandler.getStructFileB() == null);
    assert(argHandler.getMaxThreads() == -1);
    assert(argHandler.getVerbosity() == 1);
    assert(argHandler.isToUpper() == false);
    assert(argHandler.isJustDoConvert() == false);
    assert(argHandler.isJustDoSubOpt() == false);
    assert(argHandler.isJustTree() == false);
    opal.IO.Configuration[] configs = argHandler.getAdvisingConfigs();
    assert(configs.length == 1);
    System.out.println(configs[0]);
    assert(configs[0].toString().equals("DNA.260.260.69.69"));
  }

  static void test_argHandler_defaults_protein(){
    String[] argv = {"prot_2seq.fa"};
    ArgumentHandler argHandler = new ArgumentHandler(argv);
    assert(argHandler.getFileA() != null);
    assert(argHandler.getFileB() == null);
    assert(argHandler.getStructFileA() == null);
    assert(argHandler.getStructFileB() == null);
    assert(argHandler.getMaxThreads() == -1);
    assert(argHandler.getVerbosity() == 1);
    assert(argHandler.isToUpper() == false);
    assert(argHandler.isJustDoConvert() == false);
    assert(argHandler.isJustDoSubOpt() == false);
    assert(argHandler.isJustTree() == false);
    opal.IO.Configuration[] configs = argHandler.getAdvisingConfigs();
    assert(configs.length == 1);
    System.out.println(configs[0]);
    assert(configs[0].toString().equals("VTML200.45.11.42.40"));
  }

  static void test_argHandler_defaults_protein_2files(){
    String[] argv = {"prot_1seq_a.fa","prot_1seq_b.fa"};
    ArgumentHandler argHandler = new ArgumentHandler(argv);
    assert(argHandler.getFileA() != null);
    assert(argHandler.getFileB() == null);
    assert(argHandler.getStructFileA() == null);
    assert(argHandler.getStructFileB() == null);
    assert(argHandler.getMaxThreads() == -1);
    assert(argHandler.getVerbosity() == 1);
    assert(argHandler.isToUpper() == false);
    assert(argHandler.isJustDoConvert() == false);
    assert(argHandler.isJustDoSubOpt() == false);
    assert(argHandler.isJustTree() == false);
    opal.IO.Configuration[] configs = argHandler.getAdvisingConfigs();
    assert(configs.length == 1);
    System.out.println(configs[0]);
    assert(configs[0].toString().equals("VTML200.45.11.42.40"));
  }
}