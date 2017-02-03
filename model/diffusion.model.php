<?php

/**
 *
 * Diffusion Model Class
 *
 *  C/Cs = 1 - erf(x/(2*sqrt(Dt)))
 *    where C  = Concenteration
 *          Cs = Surface Level
 *          x  = x in m
 *          t  = time in seconds
 *          D  = Diffusion Constatnt (m2/s)
 *
 * @Assumptions
 *    1. Mesh centered grid is chosen
 *    2. Used central-difference form of derivatives
 *    3. Iterative Solver is used
 *
 * @author Tarun K Singhal tarun.singhal@mail.utoronto.ca
 */

use MathPHP\LinearAlgebra\Matrix;
use MathPHP\LinearAlgebra\MatrixFactory;
use MathPHP\LinearAlgebra\Vector;
use MathPHP\Functions\Special;

class DiffusionModel extends Common
{

    private $calculated_results;
    private $actual_results;
    public $error_results;
    private $STEP_SIZE;
    private $TIME_STEP;
    private $a_val = []; // a1, a2, a3
    private $b_val = []; // b1, b2, b3

    public function __construct($STEP_SIZE = STEP_SIZE, $TIME_STEP = TIME_STEP)
    {
        //set values
        $this->TIME_STEP = $TIME_STEP;
        $this->STEP_SIZE = $STEP_SIZE;

        // a values
        $this->a_val[1] = (-1.0 * THETA * DIFFUSION_CONSTANT) / ($STEP_SIZE * $STEP_SIZE);
        $this->a_val[2] = (1.0 / $TIME_STEP) + (2.0 * THETA * DIFFUSION_CONSTANT / ($STEP_SIZE * $STEP_SIZE));
        $this->a_val[3] = (-1.0 * THETA * DIFFUSION_CONSTANT) / ($STEP_SIZE * $STEP_SIZE);

        // b values
        $this->b_val[1] = (1.0 - THETA) * DIFFUSION_CONSTANT / ($STEP_SIZE * $STEP_SIZE);
        $this->b_val[2] = (1.0 / $TIME_STEP) - (2.0 * (1.0 - THETA) * DIFFUSION_CONSTANT / ($STEP_SIZE * $STEP_SIZE));
        $this->b_val[3] = (1.0 - THETA) * DIFFUSION_CONSTANT / ($STEP_SIZE * $STEP_SIZE);
    }

    public function run()
    {
        //round off mapping function
        $mat_rnd = function($x){return round($x , 5);};

        // get actual results using formula
        $this->actual_results = $this->get_actual_results();

        // get calculated results
        $this->calculated_results = $this->get_calculated_results();

        // calculate RMS error
        $this->error_results = $this->calculate_error();

        // debug
        // $this->dump($this->actual_results->map($mat_rnd)->getColumn($this->actual_results->getM() - 1));
        // $this->show('<hr>');
        // $this->dump($this->calculated_results->map($mat_rnd)->getColumn($this->actual_results->getM() - 1));
        // $this->show('<hr>');
        // $this->dump('Error:' . $this->error_results . '%');
        // die;
    }

    private function calculate_error()
    {
        //define var
        $sum = 0;
        $mat_rnd = function($x){return round($x , 5);};
        $i = intval($this->actual_results->getM() - 1);
        $j = intval($this->actual_results->getN() - 1);

        //get data
        $old = new Vector($this->actual_results->map($mat_rnd)->getColumn($j));
        $new = new Vector($this->calculated_results->map($mat_rnd)->getColumn($j));
        $diff = $old->subtract($new);
        // $this->dump($old->getVector());
        // $this->dump($new->getVector());
        // $this->dump($diff->getVector());

        //calculate absolue error
        for ($i=0; $i < $diff->getN(); $i++) {
            if($old->get($i) != 0){
                $sum += abs($diff->get($i)) / $old->get($i);
            }
        }

        //return avg percent
        return round(100 * ($sum / $diff->getN()) ,2);
    }

    private function get_calculated_results()
    {
        // total time steps deltat
        $t_steps = [];
        for ($i = 1 * $this->TIME_STEP; $i <= OVERALL_TIME; $i = $i + $this->TIME_STEP) {
            $t_steps[] = $i;
        }

        // total length part deltax
        $l_steps = [];
        for ($i = 1; $i <= (TOTAL_LENGTH / $this->STEP_SIZE); $i++) {
            $l_steps[] = $i * $this->STEP_SIZE;
        }

        // constant for a system
        $lhs_matrix = $this->initialize_matrix(count($l_steps), count($l_steps));
        foreach ($l_steps as $row => $row_d) {
            foreach ($l_steps as $col => $col_d) {
                if($row == 0){
                    if($col == 0)$lhs_matrix[$row][$col] = $this->a_val[2];
                    if($col == 1)$lhs_matrix[$row][$col] = $this->a_val[3];
                }else if($row == (count($l_steps)-1)){
                    if($col == (count($l_steps)-2))$lhs_matrix[$row][$col] = $this->a_val[1];
                    if($col == (count($l_steps)-1))$lhs_matrix[$row][$col] = $this->a_val[2] + $this->a_val[3];
                }else{
                    if(($col) == ($row - 1) && $col < (count($l_steps) - 2)){
                        $lhs_matrix[$row][$col + 0] = $this->a_val[1];
                        $lhs_matrix[$row][$col + 1] = $this->a_val[2];
                        $lhs_matrix[$row][$col + 2] = $this->a_val[3];
                    }
                }
            }
        }
        $lhs_matrix = MatrixFactory::create($lhs_matrix);

        //initial concenteration matrix for all time steps
        $conc_matrix = $this->initialize_matrix((count($l_steps) + 1), count($t_steps) + 1);
        foreach ($conc_matrix as $key => $value) {
            $conc_matrix[$key][0] = 1;
        }

        //loop for each time step
        foreach ($t_steps as $row => $time) {

            //for rhs vector for this time step
            $rhs_vector = $this->form_rhs_vector($l_steps, $conc_matrix, $row);

            //store the results
            $conc_matrix[$row + 1] = array_merge([1], $lhs_matrix->solve($rhs_vector)->getVector());
        }
        // $this->show(MatrixFactory::create($conc_matrix)->map(function($x){return round($x , 3);}));

        // return after removing initial column
        return (MatrixFactory::create($conc_matrix)->transpose());
    }

    private function form_rhs_vector($l_steps, $conc_matrix, $time_step)
    {
        /**
         * @TODO resolve issues
         */

        //for rhs vector for this time step
        $rhs_vector = array_fill(0, count($l_steps), 0.0);
        foreach ($rhs_vector as $row => $step) {
            if($row == 0){
              $rhs_vector[$row] = ($this->b_val[1] - $this->a_val[1]) * $conc_matrix[$time_step][0];
              $rhs_vector[$row] += $this->b_val[2] * $conc_matrix[$time_step][1];
              $rhs_vector[$row] += $this->b_val[3] * $conc_matrix[$time_step][2];
            }else if($row == count($l_steps)-1){
              $rhs_vector[$row] = ($this->b_val[1] * $conc_matrix[$time_step][count($l_steps)-1]);
              $rhs_vector[$row] += ($this->b_val[2] + $this->b_val[3]) * $conc_matrix[$time_step][count($l_steps)];
            }else{
              $rhs_vector[$row] = $this->b_val[1] * $conc_matrix[$time_step][$row + 0];
              $rhs_vector[$row] += $this->b_val[2] * $conc_matrix[$time_step][$row + 1];
              $rhs_vector[$row] += $this->b_val[3] * $conc_matrix[$time_step][$row + 2];
            }
        }
        //return
        return new Vector($rhs_vector);
    }

    private function get_actual_results()
    {
        // total time steps deltat
        $t_steps = [];
        for ($i = 0 * $this->TIME_STEP; $i <= OVERALL_TIME; $i = $i + $this->TIME_STEP) {
            $t_steps[] = $i;
        }

        // total length part deltax
        $l_steps = [];
        for ($i=0; $i <= (TOTAL_LENGTH / $this->STEP_SIZE); $i++) {
            $l_steps[] = $i * $this->STEP_SIZE;
        }

        // get answers
        $res = $this->initialize_matrix(count($t_steps), count($l_steps));
        foreach ($t_steps as $k => $time) {
            foreach ($l_steps as $kk => $length) {
                if($time > 0)$res[$kk][$k] = (1 - Special::errorFunction($length / (2 * sqrt($time * DIFFUSION_CONSTANT))));
            }
        }

        // form martix
        $res[0][0] = 1;
        $res = MatrixFactory::create($res);

        // return
        return $res;
    }
}

?>
