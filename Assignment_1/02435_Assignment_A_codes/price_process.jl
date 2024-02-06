using Distributions

function sample_next(previous_point)

    sample = previous_point + rand(Gamma(1.0, 2.0))*rand(Normal((5 - previous_point)*0.3, 1))
    
    if sample < 0
        sample = 0
    end

    if sample > 10
        sample = 10
    end

    return sample

end