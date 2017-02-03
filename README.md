# Particle Diffusion in 1 Dimension

A small program to analyse temporal diffusion of particles in 1-dimensional solids though FEA and comparing with actual results.

## Dependency

This project requires PHP 7+. You can install PHP 7 from `http://php.net/manual/en/install.php`


## Installation

Please follow the following steps in order to start the project.

* Create a clone wherever you want.

>$ git clone https://github.com/tarunsinghal92/particle-difussion-1d

* Install the required dependencies.

>$ composer install

* Modify the values of young modulus, etc in config.inc file which is present in includes folder.

```php
/**
 * time step (seconds)
 * sample value : 1 day = 1 * 8.64 * 10^4 seconds
 */
define('TIME_STEP', (1 * 86400));

/**
 * initial Concenteration (in percent)
 */
define('INITIAL_CONC', 1.0);

/**
 * overall time (seconds)
 * sample value : 100 day
 */
define('OVERALL_TIME', (100 * 86400));

/**
 * deltax of material (m)
 */
define('STEP_SIZE', 0.02);

/**
 * total depth 100m or 0.1m
 */
define('TOTAL_LENGTH', 0.1);

/**
 * diffusion constant (m2/s)
 */
define('DIFFUSION_CONSTANT', (0.000000000006));

/**
 * theta (Crank nicolson factor)
 */
define('THETA', (0.5));
```

## Running the Program

Run following command in `terminal`.

>$ ./run.sh

Then navigate to `http://localhost:8100` to view the results.

## License

This project is licensed under `MIT` License. See LICENSE for full license text.
