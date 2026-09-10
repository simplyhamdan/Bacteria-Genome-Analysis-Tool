from Bio import SeqIO
import matplotlib.pyplot as plt
from pathlib import Path

def visualization_menu():
    """Display the visualization menu and get the user's choice."""
    print("\nVisualization Menu:")
    print("[1] Plot Genome Sizes")
    print("[2] Plot GC Content")
    print("[3] Plot AT Content")
    print("[4] Plot Base Composition")
    print("[5] Plot K-mer Frequency")
    print("[6] Back")

    while True:
        choice = input("Enter your choice (1-6): ")

        if choice in ["1", "2", "3", "4", "5", "6"]:
            return int(choice)
        else:
            print("Invalid choice. Please enter a number between 1 and 6.")
        

def analysis_menu():
    """Display the analysis menu and get the user's choice."""
    print("\nAnalysis Menu:")
    print("[1] Individual Genome Analysis")
    print("[2] Compare Genomes")
    print("[3] K-mer Analysis")
    print("[4] Visualizations")
    print("[5] Biological Interpretation")
    print("[6] Exit")

    while True:
        choice = input("Enter your choice (1-6): ")

        if choice in ["1", "2", "3", "4", "5", "6"]:
            return int(choice)
        else:
            print("Invalid choice. Please enter a number between 1 and 6.")

def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    if len(sequence) == 0:
        return 0.0
    
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100

def load_genome_from_fasta(file_path):
    """Load a genome sequence from a FASTA file."""
    try:
        genome = SeqIO.read(file_path, "fasta")
        return genome
    except Exception as e:
        print(f"Error loading genome from {file_path}: {e}")
        return None

def load_genomes_from_folder(folder_path):
    """Load all genome sequences from FASTA files in a folder."""
    genomes = []

    data_folder = Path(folder_path)

    for file_path in data_folder.glob("*.fasta"):
        genome = load_genome_from_fasta(file_path)

        if genome is not None:
            genomes.append(genome)

    return genomes

def select_genomes(genomes):
    """Display available genomes and let the user select which ones to analyze."""

    print("\nAvailable genomes:")

    for i, genome in enumerate(genomes, start=1):
        print(f"[{i}] {genome.id}")

    while True:
        selection = input("\nSelect genomes to analyze (e.g. 1,3): ")

        selected_indices = []

        for choice in selection.split(","):
            choice = choice.strip()

            if not choice.isdigit():
                selected_indices = []
                break

            index = int(choice) - 1

            if 0 <= index < len(genomes):
                selected_indices.append(index)
            else:
                selected_indices = []
                break

        # Reject duplicate selections
        if len(selected_indices) != len(set(selected_indices)):
            print("Invalid selection. Do not select the same genome more than once.")
            continue

        if selected_indices:
            break
        else:
            print("Invalid selection. Please enter valid genome numbers.")

    selected_genomes = [genomes[i] for i in selected_indices]

    return selected_genomes

def calculate_at_content(sequence):
    """Calculate the AT content of a DNA sequence."""
    if len(sequence) == 0:
        return 0.0
    at_count = sequence.count("A") + sequence.count("T")
    return (at_count / len(sequence)) * 100

def calculate_base_composition(sequence):
    """Calculate the base composition of a DNA sequence."""
    base_counts = {
        "A": sequence.count("A"),
        "C": sequence.count("C"),
        "G": sequence.count("G"),
        "T": sequence.count("T")
    }
    return base_counts

def compare_genomes(genomes):
    """Compare the selected bacterial genomes."""

    print("\nGenome Comparison")
    print("-" * 60)

    genome_data = []

    for genome in genomes:
        sequence = str(genome.seq)

        genome_data.append({
            "id": genome.id,
            "length": len(sequence),
            "gc": calculate_gc_content(sequence),
            "at": calculate_at_content(sequence)
        })

    print(
        f"{'Genome ID':<20}"
        f"{'Length (bp)':>15}"
        f"{'GC (%)':>10}"
        f"{'AT (%)':>10}"
    )

    print("-" * 60)

    for data in genome_data:
        print(
            f"{data['id']:<20}"
            f"{data['length']:>15,}"
            f"{data['gc']:>10.2f}"
            f"{data['at']:>10.2f}"
        )

    highest_gc = max(genome_data, key=lambda x: x["gc"])
    lowest_gc = min(genome_data, key=lambda x: x["gc"])

    largest_genome = max(genome_data, key=lambda x: x["length"])
    smallest_genome = min(genome_data, key=lambda x: x["length"])

    gc_difference = highest_gc["gc"] - lowest_gc["gc"]
    size_difference = largest_genome["length"] - smallest_genome["length"]

    print("\nComparison Summary")
    print("-" * 60)

    print(
        f"Highest GC Content: {highest_gc['id']} "
        f"({highest_gc['gc']:.2f}%)"
    )

    print(
        f"Lowest GC Content: {lowest_gc['id']} "
        f"({lowest_gc['gc']:.2f}%)"
    )

    print(
        f"GC Content Difference: {gc_difference:.2f} percentage points"
    )

    print(
        f"Largest Genome: {largest_genome['id']} "
        f"({largest_genome['length']:,} bp)"
    )

    print(
        f"Smallest Genome: {smallest_genome['id']} "
        f"({smallest_genome['length']:,} bp)"
    )

    print(f"Genome Size Difference: {size_difference:,} bp")

def count_kmers(sequence, k):
    """Count k-mers in a DNA sequence."""

    if k <= 0 or k > len(sequence):
        return None

    kmer_counts = {}

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]

        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1

    return kmer_counts

def most_frequent_kmers(kmer_counts, top_n=5):
    """Return the most frequent k-mers."""

    if not kmer_counts:
        return []

    sorted_kmers = sorted(
        kmer_counts.items(),
        key=lambda item: item[1],
        reverse=True
    )

    return sorted_kmers[:top_n]

def compare_kmers(genomes, k=3, top_n=5):
    """Compare the most frequent k-mers between genomes."""
    print(f"\nMost Frequent {k}-mers Comparison")
    print("-" * 50)

    for genome in genomes:
        sequence = str(genome.seq)
        kmer_counts = count_kmers(sequence, k)
        frequent_kmers = most_frequent_kmers(kmer_counts, top_n)

        print(f"\nGenome ID: {genome.id}")
        print(f"{'K-mer':<10} {'Count':<10}")
        print("-" * 20)

        for kmer, count in frequent_kmers:
            print(f"{kmer:<10} {count:<10}")

def plot_genome_sizes(genomes):
    """Plot a bar chart of genome sizes."""
    genome_ids = [genome.id for genome in genomes]
    genome_lengths = [len(genome.seq) for genome in genomes]

    plt.figure(figsize=(10, 6))
    plt.bar(genome_ids, genome_lengths, color='skyblue')
    plt.xlabel('Genome ID')
    plt.ylabel('Genome Length (bp)')
    plt.title('Genome Sizes Comparison')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def plot_gc_content(genomes):
    """Plot a bar chart of GC content for the genomes."""
    genome_ids = [genome.id for genome in genomes]
    gc_contents = [calculate_gc_content(str(genome.seq)) for genome in genomes]

    plt.figure(figsize=(10, 6))
    plt.bar(genome_ids, gc_contents, color='lightgreen')
    plt.xlabel('Genome ID')
    plt.ylabel('GC Content (%)')
    plt.title('GC Content Comparison')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def individual_genome_analysis(genomes):
    """Display individual analysis for each selected genome."""

    for genome in genomes:
        print(f"\nGenome ID: {genome.id}")
        print(f"Genome Length: {len(genome.seq)} bp")

        # Base composition
        base_composition = calculate_base_composition(str(genome.seq))
        print("Base Composition:")

        for base, count in base_composition.items():
            print(f"  {base}: {count}")

        # GC content
        gc_content = calculate_gc_content(str(genome.seq))
        print(f"GC Content: {gc_content:.2f}%")

        # AT content
        at_content = calculate_at_content(str(genome.seq))
        print(f"AT Content: {at_content:.2f}%")

def plot_at_content(genomes):
    """Plot a bar chart of AT content for the genomes."""
    genome_ids = [genome.id for genome in genomes]
    at_contents = [calculate_at_content(str(genome.seq)) for genome in genomes]

    plt.figure(figsize=(10, 6))
    plt.bar(genome_ids, at_contents, color='lightblue')
    plt.xlabel('Genome ID')
    plt.ylabel('AT Content (%)')
    plt.title('AT Content Comparison')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def plot_base_composition(genomes):
    """Plot a bar chart of base composition for the genomes."""
    genome_ids = [genome.id for genome in genomes]
    base_compositions = [calculate_base_composition(str(genome.seq)) for genome in genomes]

    # Prepare data for plotting
    bases = ['A', 'C', 'G', 'T']
    base_counts = {base: [] for base in bases}

    for composition in base_compositions:
        for base in bases:
            base_counts[base].append(composition[base])

    # Plotting
    plt.figure(figsize=(10, 6))
    bar_width = 0.2
    index = range(len(genomes))

    for i, base in enumerate(bases):
        plt.bar([x + i * bar_width for x in index], base_counts[base], width=bar_width, label=base)

    plt.xlabel('Genome ID')
    plt.ylabel('Base Count')
    plt.title('Base Composition Comparison')
    plt.xticks(
    [x + bar_width * 1.5 for x in index],
    genome_ids,
    rotation=45
)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_kmer_frequency(genomes, k=3, top_n=5):
    """Plot the most frequent k-mers for the genomes."""
    plt.figure(figsize=(10, 6))

    for genome in genomes:
        sequence = str(genome.seq)
        kmer_counts = count_kmers(sequence, k)
        frequent_kmers = most_frequent_kmers(kmer_counts, top_n)

        kmers, counts = zip(*frequent_kmers)
        plt.bar(kmers, counts, alpha=0.5, label=genome.id)

    plt.xlabel(f'{k}-mers')
    plt.ylabel('Count')
    plt.title(f'Most Frequent {k}-mers Comparison')
    plt.legend()
    plt.tight_layout()
    plt.show()

def interpret_gc_content(genomes):
    """Interpret the GC content of the genomes."""
    print("\nGC Content Interpretation")
    print("-" * 50)

    for genome in genomes:
        gc_content = calculate_gc_content(str(genome.seq))
        interpretation = ""

        if gc_content < 40:
            interpretation = "Low GC content"
        elif 40 <= gc_content <= 60:
            interpretation = "Moderate GC content"
        else:
            interpretation = "High GC content"

        print(f"Genome ID: {genome.id} | GC Content: {gc_content:.2f}% | Interpretation: {interpretation}")

def interpret_genome_size(genomes):
    """Interpret the genome sizes of the selected genomes."""

    print("\nGenome Size Interpretation")
    print("-" * 50)

    genome_lengths = {}

    for genome in genomes:
        genome_lengths[genome.id] = len(genome.seq)

    largest_genome = max(genome_lengths, key=genome_lengths.get)
    smallest_genome = min(genome_lengths, key=genome_lengths.get)

    size_difference = (
        genome_lengths[largest_genome] - genome_lengths[smallest_genome]
    )

    print(f"Genome Size Difference: {size_difference:,} bp")
    
    print(f"Largest Genome: {largest_genome} ({genome_lengths[largest_genome]} bp)")
    print(f"Smallest Genome: {smallest_genome} ({genome_lengths[smallest_genome]} bp)")

    print("\nBiological Interpretation:")

    if size_difference < 500_000:
        print("The selected genomes have relatively similar genome sizes.")
    elif size_difference < 1_500_000:
        print("The selected genomes show a moderate difference in genome size.")
    else:
        print("The selected genomes show a substantial difference in genome size.")

    print("Genome size differences can reflect differences in genes, regulatory regions,")
    print("mobile genetic elements, and other non-coding DNA.")

    largest_length = genome_lengths[largest_genome]
    smallest_length = genome_lengths[smallest_genome]

    size_ratio = largest_length / smallest_length

    print(f"{largest_genome} has approximately {size_ratio:.2f} times the genome size of {smallest_genome}.")
    print("A larger genome does not necessarily indicate a more complex organism;")
    print("genome size can vary due to differences in gene content and genomic elements.")

def main():
    """Main function to run the bacterial genome analysis tool."""

    data_folder = Path(__file__).parent / "data"
    genomes = load_genomes_from_folder(data_folder)

    if not genomes:
        print("Failed to load genomes.")
        return

    selected_genomes = select_genomes(genomes)

    while True:
        choice = analysis_menu()

        if choice == 1:
            individual_genome_analysis(selected_genomes)

        elif choice == 2:
            compare_genomes(selected_genomes)

        elif choice == 3:
            compare_kmers(selected_genomes)

        elif choice == 4:
            while True:
                viz_choice = visualization_menu()

                if viz_choice == 1:
                    plot_genome_sizes(selected_genomes)

                elif viz_choice == 2:
                    plot_gc_content(selected_genomes)

                elif viz_choice == 3:
                    plot_at_content(selected_genomes)

                elif viz_choice == 4:
                    plot_base_composition(selected_genomes)

                elif viz_choice == 5:
                    plot_kmer_frequency(selected_genomes)

                elif viz_choice == 6:
                    break

        elif choice == 5:
            interpret_gc_content(selected_genomes)
            interpret_genome_size(selected_genomes)

        elif choice == 6:
            print("Exiting...")
            break


if __name__ == "__main__":
    main()