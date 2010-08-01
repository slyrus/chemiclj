
;;; (require 'smiles-test)
;;; (in-ns 'smiles-test)

(ns smiles-test
  (:use [chemiclj.core]
        [chemiclj.element]
        [chemiclj.smiles]))

(map #(map name %)
     (map atoms
          (bonds
           (:molecule (read-smiles-string "c=1=ccc1")))))

(:molecule (read-smiles-string "C"))
(:molecule (read-smiles-string "c"))
(:molecule (read-smiles-string "CC"))

(:molecule (read-smiles-string "C=C"))
(:molecule (read-smiles-string "C1=CC1"))

(:molecule (read-smiles-string "CC(C)CC"))

(:molecule (read-smiles-string "C1CCC2CNCCC2C1"))

(:molecule (read-smiles-string "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3"))

(:molecule (read-smiles-string "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"))

(:molecule (read-smiles-string "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O"))


(map atoms
     (bonds
       (:molecule
        (read-smiles-string
         "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3"))))

(map (fn [x] (map :_name x))
     (map atoms
          (bonds
           (:molecule
            (read-smiles-string
             "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3")))))

(def q (map atoms
          (bonds
           (:molecule
            (read-smiles-string
             "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3")))))

(map :_name (first q))

(map name (first q))

(map (fn [x] (map :_name x)) q)

(map (fn [x] (map name x)) q)

(def tam (:molecule
          (read-smiles-string
           "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3")))

(def r (:molecule (read-smiles-string
                   "C1CC1")))
