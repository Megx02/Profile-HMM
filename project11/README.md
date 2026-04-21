# Introduction
This project implements a Profile Hidden Markov Model (HMM) designed to perform sequence alignment and scoring for a protein family. The model uses a "Plan7" architecture, which incorporates Match, Insertion, and Deletion states to capture the conserved and variable regions of Multiple Sequence Alignments (MSA).


# Pseudocode

```
class PHMM:
    def read_fasta(self, filename):
        # Initialize an empty list to store sequences as strings
        sequences = []
        # Open file, parse lines, append sequences to list
        return sequences

    def column_assignments(self, sequences):
        # 'col_types' will store whether each position is 'match' or 'insertion'
        col_types = []
        num_seqs = len(sequences)
        seq_len = len(sequences[0])

        # Iterate through every column position in the alignment
        for i in range(seq_len):
            gaps = 0
            # Count gaps across all sequences for this specific column 'i'
            for j in range(num_seqs):
                if sequences[j][i] == '-':
                    gaps += 1
            
            # If >= 50% of the sequences have a gap, it is an insertion column.
            # Otherwise, it is a conserved 'match' column.
            gap_fraction = gaps / num_seqs
            if gap_fraction > 0.5:
                col_types.append("insertion")
            else:
                col_types.append("match")

        # 'L' represents the total number of Match states in our model
        L = col_types.count("match")
        return col_types, L

    def sequence_states(self, sequences, col_types, L):
        seqs_states = []
        hidden_states = set() # Use a set to avoid duplicates automatically

        for sequence in sequences:
            # Start state 
            sequence_states = ["M_0"]
            counter = 1
            
            # Map every character in the sequence to a state
            for j, char in enumerate(sequence):
                # Path A: Conserved position (Match)
                if char != '-' and col_types[j] == "match":
                    sequence_states.append(f"M_{counter}")
                    counter += 1
                
                # Path B: Missing residue in a conserved position (Deletion)
                elif char == '-' and col_types[j] == "match":
                    sequence_states.append(f"D_{counter}")
                    counter += 1
                
                # Path C: Extra residue in a variable position (Insertion)
                elif char != '-' and col_types[j] == "insertion":
                    # Insertions usually link to the previous match index
                    sequence_states.append(f"I_{counter-1}")
                
                # Path D: Gap in insertion column (Ignore, do nothing)
                elif char == '-' and col_types[j] == "insertion":
                    continue

            # Terminate the path at the End state
            sequence_states.append("“M_{L+1}")
            seqs_states.append(sequence_states)

            # Identify all unique states present across the entire alignment.
            # 'hidden_states' list serves as the index for our probability matrices.
            For seq in seq_states:
		        For state in seq:
		            # Avoid duplicate entries
			        If state not in hidden_states:
                        hidden_states.append(state)

        return seqs_states, hidden_states

	def valid_transitions(self, L):
		# Initalize a dictionary for all the possible transitions from each state
		# Since every state cannot transition to every other state
		# For example, M_i can only transition to M_(i+1), I_i, and D_(i+1)
		valid_transitions = defaultdict(list)

		# Iterate through every position from 0 to L
		for i in range(0, L+1):

			# If we haven't reached the end (i=L), there are three possible transitions
				# move forward in the alignment (match state, we move to the next position)
				# insert extra residues (insert state, we stay at the same position)
				# gap (deletion state, we move to the next position)
			if i < L:
				next_match = f"M_{i+1}"
				next_delete = f"M_{i+1}"
				insert_state = f"I_{i}"

			valid_transitions[f"M_{i}"] = [next_match, insert_state, next_delete]
			valid_transitions[f"I_{i}"] = [next_match, insert_state, next_delete]
			valid_transitions[f"D_{i}"] = [next_match, insert_state, next_delete]

			# For final position (i==L), we cannot move forward so all states transition to the end state (M_(L+1))
			else:
				valid_transitions[f"M_{i}"] = [f"M_{L+1}"]
				valid_transitions[f"I_{i}"] = [f"M_{L+1}"]
				valid_transitions[f"D_{i}"] = [f"M_{L+1}"]

	def transition_probabilities(self, seqs_states, hidded_states, valid_transitions):

```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
**Fardina Tabassum** - This project was one of the most challenging concepts I had to learn. I had to rely heavily on the materials he provided especially the HTML file to grasp what was happening in the algorithm. A major initial difficulty for me was understanding Marcus's template code that he had provided and trying to understand what he did first and how we would be inheriting the classes he provided as our code integration relied on it. There were a lot of materials we had to use simultaneously, to make sure we were on the right track. There was also some contradicting information between the HTML resources and the assignment instructions which we had to figure out which made the project challenging. Our main goal as a team was to focus on the pseudocode and then try to implement whatever we can within our initial deadline. Breaking down the pseudocode was a rigorous and time-consuming task this time around. We got together as group and spent a substantial amount of time really breaking it down piece by piece, questioning why we were doing each step and I think that personally really helped my understanding. We also had less time for this project than previous ones so we could not implement much of the Python script or notebook in time, but I think we were able to make a very comprehensive and thorough pseudocode which we can go off on to complete our implementation. 

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus# Introduction
Description of the project
