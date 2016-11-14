(ns elc5351.gerchberg_saxton
  (:use clojure.core
        [elc5351 utils fft image])
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.operators :as op]
            [clojure.core.matrix.linear    :refer [norm]])
  (:import [mikera.matrixx Matrix Matrixx]
           [mikera.vectorz AVector Vector Vectorz]
           [org.jtransforms.fft DoubleFFT_1D DoubleFFT_2D]))

(set! *warn-on-reflection* true)

(m/set-current-implementation :vectorz)

(def ^Matrix target (imread_gray "resources/Lenna.png"))
(def ^Matrix source (let [[rows cols] (m/shape target)
                          m ^Matrix (m/mul (Matrix/createRandom rows cols) 255.0)]
                      m))

(def ^Matrix source2 (let [fft_mat ^Matrix (fft2 target)]
                       (complexToReal fft_mat :type :mag)))

(defn ^Matrix gerchberg-saxton [^Matrix source ^Matrix target]
  (let []
    (loop [i ^Long (long 0)
           A ^Matrix (ifft2 (realToComplex target))]
      (if (< i 50)
        (let [B ^Matrix (composeComplex source (complexToReal A :type :phase))
              C ^Matrix (fft2 B :complex true)
              D ^Matrix (composeComplex target (complexToReal C :type :phase))
              A ^Matrix (ifft2 D)]
          (recur (unchecked-inc i) A))
        (complexToReal A :type :phase)))))
