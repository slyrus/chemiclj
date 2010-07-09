
(defproject chemiclj "0.0.0"
  :description "A computational chemistry library in clojure"
  :dependencies [[org.clojure/clojure "1.2.0-master-SNAPSHOT"]
                 [org.clojure/clojure-contrib "1.2.0-SNAPSHOT"]
                 [org.freehep/freehep-graphics2d "2.1.1"]
                 [org.freehep/freehep-graphicsio-pdf "2.1.1"]
                 [org.freehep/freehep-graphicsio-svg "2.1.1"]
                 [cdk "1.3.5"]
                 [jchempaint "3.1.2"]
                 [hiccup "0.2.6"]]
  :dev-dependencies [[swank-clojure "1.3.0-SNAPSHOT"]]
  :repositories {"freehep" "http://java.freehep.org/maven2"})

