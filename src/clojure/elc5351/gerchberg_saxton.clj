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

(def ^Matrix source (imread_gray "resources/clinton_gray.png"))
(def ^Matrix target (imread_gray "resources/trump_gray.png"))
;; (def ^Matrix source (let [[rows cols] (m/shape target)
;;                           m ^Matrix (m/mul (Matrix/createRandom rows cols) 255.0)]
;;                       m))

;; (def ^Matrix source2 (let [fft_mat ^Matrix (fft2 target)]
;;                        (complexToReal fft_mat :type :mag)))

(defn ^Matrix gerchberg-saxton [^Matrix source ^Matrix target & {:keys [iterations display]
                                                                 :or {iterations 50}}]
  (let [win_source (when display
                     (imshow source :scale 255.0))
        win_target (when display
                     (imshow target :scale 255.0))]
   (loop [i ^Long (long 0)
          A ^Matrix (fft2 source :complex false)
          result nil]
     (if (< i iterations)
       (let [ft_x_I ^Matrix (complexToReal A :type :mag)
             ft_x_P ^Matrix (complexToReal A :type :phase)
             X_next ^Matrix (composeComplex target ft_x_P)
             ift_X_next ^Matrix (ifft2 X_next :complex true)
             ift_X_next_P ^Matrix (complexToReal ift_X_next :type :phase)
             x_next ^Matrix (composeComplex source ift_X_next_P)]
         ;; (when (> i 0) (println (format "Iteration %d:" i) (mat_norm (m/sub ft_x_I result))))
         (print ". ") (flush)
         (when win_source
           (imupdate win_source (complexToReal ift_X_next :type :mag) :scale 255.0 :title (format "iteration %d" i)))
         (when win_target
           (imupdate win_target ft_x_I :scale 255.0 :title (format "iteration %d" i)))
         (recur (unchecked-inc i) (fft2 x_next :complex true) ift_X_next_P))
       (imshow result :title "result")))))


