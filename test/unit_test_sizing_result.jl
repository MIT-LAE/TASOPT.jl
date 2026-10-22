@testset "Sizing result" begin
    @testset "status contract invariants" begin
        engine_failure = [(index=ipstatic, label=:static)]

        # Exhaust the discrete status space: success must be equivalent to
        # every required solver and validity condition succeeding.
        for weight_converged in (false, true),
            balance_converged in (false, true),
            has_engine_failure in (false, true),
            takeoff_converged in (nothing, false, true),
            finite_values in (false, true)

            failures = has_engine_failure ? engine_failure : empty(engine_failure)
            result = SizingResult(weight_converged, 7, balance_converged,
                failures, takeoff_converged, finite_values)
            expected_reasons = Symbol[]
            weight_converged || push!(expected_reasons, :weight_loop_not_converged)
            balance_converged || push!(expected_reasons, :balance_not_converged)
            has_engine_failure && push!(expected_reasons, :engine_point_failure)
            takeoff_converged === false && push!(expected_reasons, :takeoff_not_converged)
            finite_values || push!(expected_reasons, :nonfinite_values)

            @test result.success == isempty(expected_reasons)
            @test result.reasons == expected_reasons
            @test result.weight_iterations == 7
            @test result.engine_failures == failures
        end
    end

    @testset "aircraft sizing integration" begin
        ac = load_default_model()
        result = size_aircraft!(ac; printiter=false)

        @test result isa SizingResult
        @test result.success
        @test result.weight_converged
        @test result.weight_iterations > 0
        @test result.balance_converged
        @test isempty(result.engine_failures)
        @test result.takeoff_converged === true
        @test result.finite_values
        @test isempty(result.reasons)
        @test ac.is_sized[1]

        failed_ac = deepcopy(ac)
        failed_ac.is_sized .= true
        failed_result = size_aircraft!(failed_ac; iter=1, initwgt=true,
            printiter=false)

        @test !failed_result.weight_converged
        @test failed_result.weight_iterations == 1
        @test :weight_loop_not_converged in failed_result.reasons
        @test !failed_result.success
        @test !failed_ac.is_sized[1]
        @test failed_ac.parg[igWMTO] > 0.0
    end
end

@testset "Horizontal-tail Newton safeguards" begin
    nominal = TASOPT._horizontal_tail_newton_step(
        [2.0 0.0; 0.0 4.0], [-2.0, 4.0], 20.0, 15.0;
        area_bounds=(1.0, 100.0), x_bounds=(0.0, 30.0),
        max_area_step=5.0, max_x_step=2.0,
    )
    @test nominal.accepted
    @test !nominal.limited
    @test nominal.dSh ≈ 1.0
    @test nominal.dxw ≈ -1.0

    ill_conditioned = TASOPT._horizontal_tail_newton_step(
        [1.0 1.0; 1.0 1.0 + eps()], [1.0, 1.0], 20.0, 15.0;
        area_bounds=(1.0, 100.0), x_bounds=(0.0, 30.0),
        max_area_step=5.0, max_x_step=2.0,
    )
    @test !ill_conditioned.accepted
    @test ill_conditioned.reason == :ill_conditioned_system
    @test isfinite(ill_conditioned.condition_metric)
    @test all(isfinite, (ill_conditioned.dSh, ill_conditioned.dxw))

    bounded = TASOPT._horizontal_tail_newton_step(
        [1.0 0.0; 0.0 1.0], [100.0, -100.0], 20.0, 15.0;
        area_bounds=(1.0, 100.0), x_bounds=(0.0, 30.0),
        max_area_step=5.0, max_x_step=2.0,
    )
    @test bounded.accepted
    @test bounded.limited
    @test abs(bounded.dSh) <= 5.0
    @test abs(bounded.dxw) <= 2.0
    @test 1.0 <= bounded.Sh <= 100.0
    @test 0.0 <= bounded.xwbox <= 30.0
end


@testset "Weight-fraction validity" begin
    @test TASOPT._weight_fraction_failure(prevfloat(1.0), :test) === nothing
    for value in (1.0, 1.01, Inf, NaN)
        failure = TASOPT._weight_fraction_failure(value, :test)
        @test failure isa WeightFractionFailure
        @test failure.classification == (isfinite(value) ? :at_or_above_one : :nonfinite)
        result = SizingResult(true, 1, true, NamedTuple[], true, true, failure)
        @test !result.success
        @test :invalid_weight_fraction_sum in result.reasons
    end
end
