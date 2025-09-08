1. What is Double-Double Precision?

A double-double number is a software-emulated floating-point type that extends the precision of a standard IEEE-754 double (64-bit, ~53 bits of mantissa ≈ 16 decimal digits) by representing one logical number as the unevaluated sum of two doubles:

x = x_\text{hi} + x_\text{lo}, \quad |x_\text{lo}| \leq \tfrac{1}{2} \, \text{ulp}(x_\text{hi})
	•	hi holds the main value (nearest double).
	•	lo carries the round-off error or residual.
	•	Together, they give ~106 bits of mantissa (~31–32 decimal digits).

This is not a new IEEE format; it’s a trick relying on careful use of standard double arithmetic.

⸻

2. Core Tools: Error-Free Transforms

The whole scheme relies on algorithms that decompose floating-point results into exact parts.
	•	TwoSum(a, b) → (s, e)
Computes sum s = a + b and error e such that
a + b = s + e \quad \text{exactly.}
	•	TwoProd(a, b) → (p, e)
Computes product p = a·b and error e such that
a \cdot b = p + e \quad \text{exactly.}
(relies on fused multiply-add (FMA) when available for speed/accuracy).

These guarantee error-free transformations within IEEE double, critical for DD arithmetic.

⸻

3. Double-Double Addition

Given two DD numbers x = (x_h, x_l), y = (y_h, y_l):
	1.	(s, e) = TwoSum(x_h, y_h)
— Add the high parts, capture residual.
	2.	t = x_l + y_l
— Add the low parts.
	3.	(s2, e2) = TwoSum(s, t)
— Fold in low parts.
	4.	Result:
z_h = s2, \quad z_l = e + e2

After a possible renormalization step (making sure |z_l| < ulp(z_h)), this gives ~106 bits of precision.

⸻

4. Double-Double Multiplication

Given two DD numbers x = (x_h, x_l), y = (y_h, y_l):
	1.	(p, e) = TwoProd(x_h, y_h)
— Main product and its rounding error.
	2.	p1 = p
— Start with main product.
	3.	p2 = (x_h * y_l) + (x_l * y_h)
— Cross terms (each smaller).
	4.	Sum them carefully:
(s, e2) = TwoSum(p1, p2)
	5.	Accumulate all errors:
lo = e + e2 + (x_l * y_l)
	6.	Renormalize:
(hi, lo2) = TwoSum(s, lo)

Final result:
z_h = hi, \quad z_l = lo2

This achieves double-double precision multiplication with error bounded to about one ulp of the 106-bit result.

⸻

5. Why It Works
	•	IEEE-754 double has guard digits and FMA that let you recover round-off exactly.
	•	By keeping the small remainder (lo) separate, you retain ~53 + 53 = 106 bits.
	•	DD math is significantly slower (≈ 10–20× a plain double multiply), but much faster than arbitrary-precision libraries.

⸻

6. Other Operations
	•	Division: Compute hi = x_h / y_h, refine via Newton iteration with DD multiply.
	•	Square root: Similar, use initial double guess, refine with Newton.
	•	Normalization: After every op, adjust so that hi is the rounded double and lo is the small correction.

⸻

7. Applications
	•	Fractal zooming (Mandelbrot/Julia)
	•	High-precision geometry kernels
	•	Numerical analysis (compensated algorithms)
	•	Validating IEEE quad-precision implementations
