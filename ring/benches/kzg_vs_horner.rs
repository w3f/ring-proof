use ark_ec::Group;
use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_std::{test_rng, UniformRand};
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};


fn horner<F: Field>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let mut group = c.benchmark_group("horner");
    for log_d in 10..=12 {
        let d = 2usize.pow(log_d);
        group.bench_with_input(BenchmarkId::from_parameter(d), &d, |b, &d| {
            let poly = DensePolynomial::rand(d, rng);
            b.iter_with_setup(
                || F::rand(rng),
                |z| poly.evaluate(&z),
            )
        });
    }
    group.finish();
}

fn kzg<G: Group>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let mut group = c.benchmark_group("kzg");
    group.bench_function("kzg", |b|
        b.iter_with_setup(
            || (G::ScalarField::from(u128::rand(rng)), G::rand(rng)),
            |(x, g)| g * x,
        ),
    );
    group.finish();
}

criterion_group!(benches, horner::<ark_bls12_381::Fr>, kzg::<ark_bls12_381::G1Projective>);
criterion_main!(benches);