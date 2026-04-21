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






	# Function that finds the emission probabilities of each state in the model
	def emission_probabilites(col_types, seqs, emissions, L, b):

	  # Initialize emission probability dictionary
	  emission_probs = defaultdict(dict)
	     
	  # Get emission counts for all match states
	  idx = 1
	  # Looping over col_types ensures we reach every aa in a sequence
	  for state in enumerate(col_types):
	    for seq in seqs:
	    # if we are at a match state add counts for each amino acid encountered in the seqs at that state
	     if state == "match":
	      emission_probs[f"M_{idx}"][seq[i]] += 1
	      # Add 1 to the index to keep total match state (L) integrity
	
	  # Convert emissions to probabilities
	  for state in emission_probs.keys():
	    s = sum(emission_probs[state].values())
	    for emission in emission_probs[state].keys():
	      # Update each emission in each state with count / sum of all counts in that state (including pseudocounts to prevent division with 0)
	      emission_probs[state][emission] = emission_probs[state][emission] + 1b / s + 20b
	  
	  # Now, get emission counts for all insertion states, which is just the background frequency of amino acids across all sequences
	  background_probs = defaultdict()
	    for seq in seqs:
	      for emission in seq:
	        if emission != '-':
	          background_probs[emission] += 1
	  
	  # Convert the counts we have now to probabilities
	  total = sum(background_probs.values())
	  for emission in background_probs.keys():
	    background_probs[emission] = background_probs[emission] + 1b / total + 20b
	    
	  # Assign background probs to all insertion states
	  for i in range(self.L + 1):
	    emission_probs[f"I_{i}"] = background_probs.copy()
	  
	  return emission_probs
 
```

# Successes
The biggest success of this project was our groups ability to collaborate to come up with a roadmap for this project. We don't have full implementation complete, but that's because we spent so much time making sure we all understood how to create the pHMM topology from MSA and had pseudocode that could be easily transferred to code. All of the pseudocode outlined above was done as a group over many hours, collecting information from different resources and problem solving in real time. There were a lot of steps and moving parts in this project and our ability to collaborate made it much easier to synthesize the material and decide on what was most important. Our group really honed in on the pHMM topology creation, and while we haven't worked through using Viterbi and Forward-Backward on a new amino acid sequence, the construction of our model as dictionaries should make the integration with the provided HMM.py module relatively smooth. Overall, the plan and pseudocode that we've come up with should be a sufficient basis for performing all the tasks of the project if used for implementation. 

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
**Fardina Tabassum** - This project was one of the most challenging concepts I had to learn. I had to rely heavily on the materials he provided especially the HTML file to grasp what was happening in the algorithm. A major initial difficulty for me was understanding Marcus's template code that he had provided and trying to understand what he did first and how we would be inheriting the classes he provided as our code integration relied on it. There were a lot of materials we had to use simultaneously, to make sure we were on the right track. There was also some contradicting information between the HTML resources and the assignment instructions which we had to figure out which made the project challenging. Our main goal as a team was to focus on the pseudocode and then try to implement whatever we can within our initial deadline. Breaking down the pseudocode was a rigorous and time-consuming task this time around. We got together as group and spent a substantial amount of time really breaking it down piece by piece, questioning why we were doing each step and I think that personally really helped my understanding. We also had less time for this project than previous ones so we could not implement much of the Python script or notebook in time, but I think we were able to make a very comprehensive and thorough pseudocode which we can go off on to complete our implementation. 

## Other member
**Connor Crawford** - The scope of this project was what was most difficult for me to deal with. There were a lot of outlined expectations along with integration with the HMM code we were provided that it was hard for to orient towards the material. Thankfully, once our group started meeting they were able to talk me through their perspectives on the project and point me towards some helpful resources. Once we started planning as a group I felt much more comfortable with the material, and more specifically, the steps for constructing the pHMM topology became much clearer. I still feel as though I'm a little confused on how the code integration and more complex analysis of this project will work but I definitely feel confident with the base steps of the project. When doing the planning/pseudocode we spent many hours meeting as a group, and went through each function together - mapping out what it needed to do, and how we thought we could go about doing it. I thought it was really beneficial to get everyone's perspective and talk through problems together, I know that I would not have been able to understand the project half as well if it weren't for this collaboration. 

# Generative AI Appendix
As per the syllabus

