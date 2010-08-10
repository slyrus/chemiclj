;;; file: test/smiles_test.clj
;;;
;;; Copyright (c) 2010 Cyrus Harmon (ch-lisp@bobobeach.com) All rights
;;; reserved.
;;;
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
;;;
;;;   * Redistributions of source code must retain the above copyright
;;;     notice, this list of conditions and the following disclaimer.
;;;
;;;   * Redistributions in binary form must reproduce the above
;;;     copyright notice, this list of conditions and the following
;;;     disclaimer in the documentation and/or other materials
;;;     provided with the distribution.
;;;
;;; THIS SOFTWARE IS PROVIDED BY THE AUTHOR 'AS IS' AND ANY EXPRESSED
;;; OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
;;; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;;; ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
;;; DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
;;; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
;;; GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;;; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
;;; WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


;;; (require 'smiles-test)
;;; (in-ns 'smiles-test)

(ns smiles-test
  (:use [chemiclj.core]
        [chemiclj.element]
        [chemiclj.smiles]))

(read-smiles-string "C")
(read-smiles-string "CC")
(read-smiles-string "C=C")
(read-smiles-string "C1=CC1")
(read-smiles-string "CC(C)CC")
(read-smiles-string "C1CCC2CNCCC2C1")
(read-smiles-string "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3")
(read-smiles-string "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O")
(read-smiles-string "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O")

(def *molecules*
     (reduce (fn [acc [name smiles]]
               (into acc {name (name-molecule (read-smiles-string smiles) name)}))
             {}
             {"methane" "C"
              "ethane" "CC"
              "propane" "CCC"
              "butane" "CCCC"
              "pentane" "CCCCC"
              "hexane" "CCCCCC"
              "l-alanine" "C[C@@H](C(=O)O)N"
              "valine" "CC(C)C(C(=O)O)N"
              "paroxetine" "C1CNCC(C1C2=CC=C(C=C2)F)COC3=CC4=C(C=C3)OCO4"
              "tamoxifen" "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3"
              "anastrozole" "CC(C)(C#N)C1=CC(=CC(=C1)CN2C=NC=N2)C(C)(C)C#N"
              "acetominophen" "CC(=O)NC1=CC=C(C=C1)O"
              "morphine" "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"
              "estradiol" "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"
              "dicyclohexyl" "C1CCCCC1C2CCCCC2"
              "spiro[5.5]undecane" "C12(CCCCC1)CCCCC2"
              "benzene" "c1ccccc1"
              "indane" "c1ccc2CCCc2c1"
              "furan" "c1occc1"
              "fluoroform" "C(F)(F)F"
              "vanillin" "O=Cc1ccc(O)c(OC)c1"
              "thiamin" "OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2"
              "oxytocin" "[H][C@]1(NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O)[C@@H](C)CC"
              "sildenafil" "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1"
              "capsaicin" "COc1cc(CNC(=O)CCCC\\C=C\\C(C)C)ccc1O"
              "tricky" "[C@@]123[C@H](C(C=C3)(C)C)CC[C@@](C1)(CCC2)C"
              "6-amino-2-ethyl-5-(aminomethyl)-1-hexanol" "OCC(CC)CCC(CN)CN"}))

(defn get-molecule [name]
  (get *molecules* name))

(read-smiles-string "C[C@](Br)(Cl)I") 
(read-smiles-string "C[C@H](Br)I")
(read-smiles-string "N[C@](Br)(O)C")
(read-smiles-string "N[C@@](Br)(C)O")

