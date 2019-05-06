class Transposon():
    """Handle a transposon object."""

    def __init__(self, seq_name, first, last, score):
        """Construct for the transposon instance."""
        self.sequence_name = seq_name
        self.first = int(first)
        self.last = int(last)
        self.score = int(score)

    def set_sequence_name(self, seq_name):
        """Set sequence name."""
        self.sequence_name = seq_name

    def set_first(self, first):
        """Set first."""
        self.first = first

    def set_last(self, last):
        """Set last."""
        self.last = last

    def set_score(self, score):
        """Set score."""
        self.score = score

    def is_overlap(self, transposon):
        """Verify if there's an overlap with another transposon."""
        if self.first <= transposon.last <= self.last:
            return True
        elif self.first <= transposon.first <= self.last:
            return True
        else:
            return False

    def get_overlap(self, transposon):
        """Get the size of overlap with another transposon."""
        if not self.is_overlap(transposon):
            return 0
        end = min(self.last, transposon.last)
        start = max(self.first, transposon.first)
        return end - start + 1

    def __len__(self):
        """Get size of the element."""
        return self.last - self.first + 1

    def __str__(self):
        """Write the transposon data as string"""
        return self.sequence_name + '\t' + self.first + '\t' + self.last
