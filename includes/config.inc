<?php

/**
 * url
 */
define('SITEURL', 'particle-difussion-1d/');

/**
 * time step (seconds)
 * sample value : 1 day = 1 * 8.64 * 10^6 seconds
 */
define('TIME_STEP', (5 * 8.64 * 1000000));

/**
 * initial Concenteration (in percent)
 */
define('INITIAL_CONC', 1.0);

/**
 * overall time (seconds)
 * sample value : 100 day
 */
define('OVERALL_TIME', 100 * 8.64 * 1000000);

/**
 * deltax of material (m)
 */
define('STEP_SIZE', 0.001);

/**
 * Number of 1D elements
 */
define('NUM_ELEMENTS', 100);

/**
 * diffusion constant (m2/s)
 */
define('DIFFUSION_CONSTANT', (0.000000000006));

/**
 * theta (Crank nicolson factor)
 */
define('THETA', (0.5));

?>
