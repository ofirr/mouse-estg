function [ newMS ] = update_microsatellite( MS, T, Nodes, parentType, rnd, rep )
% updates the microsatellite repeat number

    rand_seed = RandStream('mlfg6331_64');

    % access to simulation options
    global simul_options;

    % access to microsatellite mutation transition table
    global ms_mutation_transition_prob;

    % mapping between index and ms repeat length
    global ms_idx_rptlen_mapping;

    % check if MS has all repeat number -1, which means it requires
    % reinitialization
    if sum( MS ~= -1 ) == 0

        % access to om6 microsatellite ids and repeat numbers read
        global om6_ms_alleles;

        % reinitialize using the ms repeat numbers with biallelic
        % om6_ms_alleles{1} : initial repeat lengths for mono-allelic
        MS = om6_ms_alleles{1}(:,2)';

    end

    % don't induce mutations if the option is turned off
    if ~simul_options.addMutations
        newMS = MS;
        return;
    end

    % get the time that elapsed since the previous cell division
    parent_time = Nodes{parentType}(rnd).InternalStates.('Time');
    time_factor = T - parent_time;
    if min(MS)==5
    for i = 1:length(MS)

        try
            % convert ms repeat length to index
            idx = find(ms_idx_rptlen_mapping == MS(i));

            % using the index, get probability distribution
            probs = ms_mutation_transition_prob(idx, :);
        catch
            fprintf("error: %d\n", i);
        end

        %fixme: with rand_seed, barely changes

        % update all transition probabilities by the time factor
        if max(probs) < 1
            probs = [probs(1) probs(2:length(probs))*0.01*time_factor];
        end

        % get index that corresponds to ms repeat length using probability
        % distribution
        new_idx = datasample(...
            1:length(probs), ...
            1, ...
            'Replace', true, ...
            'Weights', probs ...
        );

        % convert index to ms repeat length
        newMS(i) = ms_idx_rptlen_mapping(new_idx);

    end
    if sum(MS~=newMS) > 1
        % roll the dice for multi site deletion
        probs = ms_mutation_transition_prob(1, :);
        probs = [probs(1) probs(2:length(probs))*100000];
        deletion_threshold_index = 2;
        % get index that corresponds to the deletion probability 1..i indicates 
        % no deletion and i..n indicates a deletion of n scars
        deleted_n = datasample(...
            1:length(probs), ...
            1, ...
            'Replace', true, ...
            'Weights', probs ...
        );
        if deleted_n > deletion_threshold_index
            if length(MS)-2-deleted_n > 0
                first_deleted_scar_index = randi(length(MS)-3-deleted_n);
                for i = first_deleted_scar_index+1:first_deleted_scar_index+deleted_n+1
                    newMS(i) = 28;  % assign deleted character
                end
            end
        end
    end
    else
        newMS = MS;
    end
end
