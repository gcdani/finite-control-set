function output = load_torque_create( input, start_time, end_time, sim_time )

    if sim_time >= start_time && sim_time <= end_time
        output = input;
    else
        output = 0;
    end


end

