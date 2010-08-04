
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
             {"valine" "CC(C)C(C(=O)O)N"
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
              "thiamin" "OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2"}))

(defn get-molecule [name]
  (get *molecules* name))

