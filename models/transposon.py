class Transposon():
    """Handle a transposon object."""

    def __init__(self, seq_name, first, last):
        """Construct for the transposon instance."""
        self.sequence_name = seq_name
        self.first = int(first)
        self.last = int(last)

    def set_sequence_name(self, seq_name):
        """Set sequence name."""
        self.sequence_name = seq_name

    def set_first(self, first):
        """Set first."""
        self.first = first

    def set_last(self, last):
        """Set last."""
        self.last = last

    def __len__(self):
        """Get size of the element."""
        return self.last - self.first + 1

    def __str__(self):
        """Write the transposon data as string"""
        return self.sequence_name + '\t' + self.first + '\t' + self.last
