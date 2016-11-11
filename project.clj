(defproject ELC5351_MS "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :source-paths ["src/clojure/ELC5351_MS/"]
  :dependencies [[org.clojure/clojure      "1.8.0"]
                 [org.clojure/tools.cli    "0.3.5"]
                 [net.mikera/core.matrix   "0.56.0"]
                 [net.mikera/vectorz       "0.64.0"]
                 [net.mikera/vectorz-clj   "0.45.0"]
                 ;; [clatrix "0.5.0"]
                 ;; [com.github.wendykierp/JTransforms "3.1"]
                 [incanter/incanter-charts "1.5.7"]
                 [incanter/incanter-core   "1.5.7"]
                 ;; [incanter/incanter-charts "1.9.0"]
                 ;; [incanter/incanter-core   "1.9.0"]
                 [org.clojure/core.async   "0.2.385"]
                 [commons-validator        "1.4.1"]
                 [org.apache.commons/commons-math3 "3.5"]
                 [org.clojure/tools.logging "0.3.1"]]
  :jvm-opts ^:replace []
  :repl-options {:port 5351})

