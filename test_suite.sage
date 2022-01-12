attach("mf_secret_sharing.sage")

def test_recon_random(level, weight, prime, num_senders, r):
    ring = GF(prime)
    M = ModularForms(level, weight, base_ring = ring)
    print("Using MFspace %s"%M)
    pub = create_public_data(M, num_senders, r)
    print("Using public data:\n%s"%pub)
    votes = [ring(randint(0,1)) for ix in range(0,num_senders)]
    sum_votes = sum(votes)
    print("Vote to be cast:\n%s"%sum_votes)

    #senders send (partial) votes to auditors
    auditors = [R.zero() for ix in range(0,M.sturm_bound()+1)]
    for sender in range(0,num_senders):
        for auditor in range(0,len(auditors)):
            #auditors collect (partial) votes
            auditors[auditor] += send_to_auditor(votes[sender], sender, auditor, pub)

    #auditors send (partial) votes to senders
    sender_recv=[[] for ix in range(0,num_senders)]
    for auditor in range(0,len(auditors)):
        for sender in range(0,num_senders):
            #need id because order is important, and comms might be async
            sender_recv[sender].append((auditor,auditors[auditor]))

    #senders reconstruct vote sum
    sum_votes_recon = [R.zero() for ix in range(0,num_senders)]
    for sender in range(0,num_senders):
        sum_votes_recon[sender] = reconstruct(sender, sender_recv[sender], pub)
        assert sum_votes_recon[sender] == sum_votes
