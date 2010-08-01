
(ns smiles-test
  (:use [chemiclj.core]
        [chemiclj.element]
        [chemiclj.smiles]))

(map #(map name %)
     (map atoms
          (bonds
           (:molecule (chemiclj.smiles/read-smiles-string "c=1=ccc1")))))

(:molecule (chemiclj.smiles/read-smiles-string "C"))
(:molecule (chemiclj.smiles/read-smiles-string "c"))
(:molecule (chemiclj.smiles/read-smiles-string "CC"))

(:molecule (chemiclj.smiles/read-smiles-string "C=C"))
(:molecule (chemiclj.smiles/read-smiles-string "C1=CC1"))

(:molecule (chemiclj.smiles/read-smiles-string "CC(C)CC"))

(:molecule (chemiclj.smiles/read-smiles-string "C1CCC2CNCCC2C1"))

 (:molecule (chemiclj.smiles/read-smiles-string "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3"))

(:molecule (chemiclj.smiles/read-smiles-string "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"))

(:molecule (chemiclj.smiles/read-smiles-string "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"))

