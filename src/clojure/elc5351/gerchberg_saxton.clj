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

(def ^Matrix source (imread_gray "resources/Lenna.png"))
(def ^Matrix target (complexToReal (fft2 source) :type :mag))
(def ^Matrix source (imread_gray "resources/intensity_512x512.png"))
;; (def ^Matrix source (let [[rows cols] (m/shape target)
;;                           m ^Matrix (m/mul (Matrix/createRandom rows cols) 255.0)]
;;                       m))

;; (def ^Matrix source2 (let [fft_mat ^Matrix (fft2 target)]
;;                        (complexToReal fft_mat :type :mag)))

(defn ^Matrix gerchberg-saxton [^Matrix source ^Matrix target & {:keys [iterations]
                                                                 :or {iterations 50}}]
  (let []
    (loop [i ^Long (long 0)
           x ^Matrix (realToComplex source)
           win (imshow (complexToReal (fft2 source) :type :phase) :scale 255.0 :title :phase)]
      (if (< i iterations)
        (let [A ^Matrix (fft2 x :complex true :copy false)
              A_phase ^Matrix (complexToReal A :type :phase)
              new_FFT ^Matrix (composeComplex target A_phase)
              C       ^Matrix (ifft2 new_FFT :complex true :copy false)
              C_phase ^Matrix ()
              D ^Matrix (composeComplex source (complexToReal C :type :phase))]
          (if (= i 0)
            (imshow C :scale 255.0 :title (format "iteration %d" i))
            (imupdate ))

          (print ". ") (flush)
          (recur (unchecked-inc i) D))
        (do
          (imshow (complexToReal x :type :mag)   :scale 255.0 :title "Using Magnitude")
          (imshow (complexToReal x :type :phase) :scale 255.0 :title "Using Phase")
          (complexToReal x :type :phase))))))

(defn ^Matrix gerchberg-saxton [^Matrix source ^Matrix target & {:keys [iterations]
                                                                 :or {iterations 50}}]
  (let [A ^Matrix (ifft2 target :complex false)
        win (imshow (complexToReal A :type :phase) :scale 255.0 :title :phase)]
    (println "size A" (m/shape A))
    (loop [i ^Long (long 0)
           A A
           result nil]
      (if (< i iterations)
        (let [B ^Matrix (composeComplex source (complexToReal A :type :phase))
              C ^Matrix (fft2  B :complex true :copy false)
              D ^Matrix (composeComplex target (complexToReal C :type :phase))
              A ^Matrix (ifft2 D :complex true :copy false)]
          (if (= i 0)
            (imshow (complexToReal C :type :mag) :scale 255.0 :title (format "iteration %d" i))
            (imupdate win (complexToReal C :type :mag) :scale 255.0 :title (format "iteration %d" i)))
          (print ". ") (flush)
          (recur (unchecked-inc i) A (complexToReal C :type :mag)))
        (do
          (imshow result   :scale 255.0 :title "Using Magnitude")
          result )))))
